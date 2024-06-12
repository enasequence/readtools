/*
 * Copyright 2010-2021 EMBL - European Bioinformatics Institute
 * Licensed under the Apache License, Version 2.0 (the "License"); you may not use this
 * file except in compliance with the License. You may obtain a copy of the License at
 * http://www.apache.org/licenses/LICENSE-2.0
 * Unless required by applicable law or agreed to in writing, software distributed under the
 * License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
 * CONDITIONS OF ANY KIND, either express or implied. See the License for the
 * specific language governing permissions and limitations under the License.
 */
package uk.ac.ebi.ena.readtools.cram.ref;

import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.cram.io.InputStreamUtils;
import htsjdk.samtools.cram.ref.CRAMReferenceSource;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.lang.ref.Reference;
import java.lang.ref.SoftReference;
import java.net.MalformedURLException;
import java.net.URL;
import java.net.URLConnection;
import java.time.Duration;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicLong;
import java.util.stream.Collectors;
import org.springframework.retry.backoff.ExponentialBackOffPolicy;
import org.springframework.retry.policy.SimpleRetryPolicy;
import org.springframework.retry.support.RetryTemplate;
import uk.ac.ebi.ena.readtools.cram.common.Utils;

/**
 * A central class for automated discovery of reference sequences. The algorithm is expected similar
 * to that of samtools:
 *
 * <ul>
 *   <li>Search in memory cache by sequence name.
 *   <li>Use local fasta file is supplied as a reference file and cache the found sequence in
 *       memory.
 *   <li>Try REF_CACHE env variable.
 *   <li>Try all entries in REF_PATH. The default value is the EBI reference service.
 *   <li>Try @SQ:UR as a URL for a fasta file with the fasta index next to it.
 * </ul>
 */
public class ENAReferenceSource implements CRAMReferenceSource {
  public interface LoggerWrapper {
    static final Log log = Log.getInstance(ENAReferenceSource.class);

    public default void debug(final Object... messageParts) {
      log.debug(messageParts);
    }

    public default void warn(final Object... messageParts) {
      log.warn(messageParts);
    }

    public default void error(final Object... messageParts) {
      log.error(messageParts);
    }

    public default void info(final Object... messageParts) {
      log.info(messageParts);
    }
  }

  private static final String JAVA_IO_TMPDIR_PROPERTY_NAME = "java.io.tmpdir";
  private static final String USER_HOME_PROPERTY_NAME = "user.home";

  private static final int REF_BASES_TO_CHECK_FOR_SANITY = 1000;

  private static String REF_CACHE = System.getenv("REF_CACHE");
  private static String XDG_CACHE = System.getenv("XDG_CACHE_HOME");
  private static String DEFAULT_CACHE = ".cache/hts-ref/%2s/%2s/%s";

  private static String REF_PATH = System.getenv("REF_PATH");

  private List<PathPattern> refPatterns = new ArrayList<PathPattern>();
  private List<PathPattern> cachePatterns = new ArrayList<PathPattern>();

  Map<String, Integer> disk_map = new HashMap<>();
  private AtomicInteger disk_counter = new AtomicInteger();
  private AtomicLong disk_sz = new AtomicLong();
  private AtomicLong disk_spent = new AtomicLong();

  Map<String, Integer> download_map = new HashMap<>();
  private AtomicInteger download_counter = new AtomicInteger();
  private AtomicLong download_sz = new AtomicLong();
  private AtomicLong download_spent = new AtomicLong();

  private AtomicInteger mem_counter = new AtomicInteger();

  int getRemoteFetchCount() {
    return download_counter.get();
  }

  int getDiskFetchCount() {
    return disk_counter.get();
  }

  int getMemFetchCount() {
    return mem_counter.get();
  }

  /*
   *
   *  Initialisation is according to Samtools' misc/seq_cache_populate.pl
   *
   */
  {
    if (REF_PATH == null) REF_PATH = "https://www.ebi.ac.uk/ena/cram/md5/%s";

    if (null != REF_CACHE) {
      cachePatterns.add(new PathPattern(REF_CACHE));
    } else if (null != XDG_CACHE) {
      cachePatterns.add(new PathPattern(XDG_CACHE + "/" + DEFAULT_CACHE));
    } else if (null != System.getProperty(USER_HOME_PROPERTY_NAME)) {
      cachePatterns.add(
          new PathPattern(System.getProperty(USER_HOME_PROPERTY_NAME) + "/" + DEFAULT_CACHE));
    } else if (null != System.getProperty(JAVA_IO_TMPDIR_PROPERTY_NAME)) {
      cachePatterns.add(
          new PathPattern(System.getProperty(JAVA_IO_TMPDIR_PROPERTY_NAME) + "/" + DEFAULT_CACHE));
    }

    refPatterns.addAll(
        splitRefPath(REF_PATH).stream().map(e -> new PathPattern(e)).collect(Collectors.toList()));
  }

  // TODO more strict regexp
  static List<String> splitRefPath(String paths) {
    String r = "(?i)(?<!(https?|ftp))(?!:([0-9]+))(:|;)";
    return Arrays.stream(paths.split(r)).collect(Collectors.toList());
  }

  public List<String> getRefPathList() {
    return Collections.unmodifiableList(
        refPatterns.stream().map(e -> e.toString()).collect(Collectors.toList()));
  }

  public List<String> getRefCacheList() {
    return Collections.unmodifiableList(
        cachePatterns.stream().map(e -> e.toString()).collect(Collectors.toList()));
  }

  private LoggerWrapper log = new LoggerWrapper() {};
  private int downloadTriesBeforeFailing = 2;

  /*
   * In-memory cache of ref bases by sequence name. Garbage collector will
   * automatically clean it if memory is low.
   */
  private Map<String, Reference<byte[]>> cacheW =
      new ConcurrentHashMap<String, Reference<byte[]>>();

  public ENAReferenceSource() {
    ;
  }

  public ENAReferenceSource(String... cachePatterns) {
    this.cachePatterns.clear();
    this.cachePatterns.addAll(
        Arrays.stream(cachePatterns).map(e -> new PathPattern(e)).collect(Collectors.toList()));
  }

  public void setLoggerWrapper(LoggerWrapper log) {
    this.log = log;
  }

  public void clearMemCache() {
    cacheW.clear();
  }

  private byte[] findInMemCache(String md5) {
    Reference<byte[]> r = cacheW.get(md5);
    if (r != null) {
      byte[] bytes = r.get();
      if (bytes != null) {
        log.debug(
            String.format(
                "% 6d Reference found in memory cache by md5: %s",
                mem_counter.incrementAndGet(), md5));
        return bytes;
      }
    }
    return null;
  }

  private void addToMemCache(String md5, byte[] bytes) {
    cacheW.put(md5, new SoftReference<byte[]>(bytes));
  }

  private byte[] findInFileCache(String md5) throws IOException {
    long start = System.currentTimeMillis();
    for (PathPattern pathPattern : cachePatterns) {
      File file = new File(pathPattern.format(md5));
      if (file.exists()) {
        byte[] data = loadFromPathWithFileCleanupRetry(file, md5);
        log.debug(
            String.format(
                ".% 5d Reference found on disk cache at the location %s sz:%d, total: %d, spent: %d, total spent: %d, attempt %d",
                disk_counter.incrementAndGet(),
                file.getPath(),
                data.length,
                disk_sz.addAndGet(data.length),
                System.currentTimeMillis() - start,
                disk_spent.addAndGet(System.currentTimeMillis() - start),
                disk_map.merge(md5, (Integer) 1, (v1, v2) -> v1 + v2)));
        return data;
      }
    }

    return null;
  }

  @Override
  public synchronized byte[] getReferenceBases(SAMSequenceRecord record, boolean tryNameVariants) {
    byte[] bases = findBases(record);
    if (bases == null) return null;

    String md5 = record.getAttribute(SAMSequenceRecord.MD5_TAG);
    if (md5 == null) {
      md5 = Utils.calculateMD5String(bases);
      record.setAttribute(SAMSequenceRecord.MD5_TAG, md5);
    }

    return bases;
  }

  @Override
  public byte[] getReferenceBasesByRegion(
      SAMSequenceRecord record, int zeroBasedStart, int requestedRegionLength) {
    byte[] bases = getReferenceBases(record, false);

    int endIndex = calculateEndIndex(zeroBasedStart, requestedRegionLength, bases.length);
    return Arrays.copyOfRange(bases, zeroBasedStart, endIndex);
  }

  private int calculateEndIndex(int start, int length, int maxLength) {
    int endIndex = start + length;
    return Math.min(endIndex, maxLength);
  }

  private byte[] readBytesFromFile(File file, int offset, int len) throws IOException {
    long size = file.length();
    if (size < offset || len < 0) {
      log.warn(
          String.format(
              "Ref request is out of range: %s, size=%d, offset=%d, len=%d",
              file.getAbsolutePath(), size, offset, len));
      return new byte[] {};
    }

    byte[] data = new byte[(int) Math.min(size - offset, len)];
    try (FileInputStream fis = new FileInputStream(file);
        DataInputStream dis = new DataInputStream(new BufferedInputStream(fis)); ) {
      dis.skip(offset);
      dis.readFully(data);
      dis.close();
      return data;
    }
  }

  public synchronized ReferenceRegion getRegion(
      SAMSequenceRecord record, int start_1based, int endInclusive_1based) throws IOException {

    String md5 = record.getAttribute(SAMSequenceRecord.MD5_TAG);
    if (md5 == null) return null;

    { // check cache by md5:
      byte[] bases = findInMemCache(md5);
      if (bases != null) {
        return ReferenceRegion.copyRegion(
            bases,
            record.getSequenceIndex(),
            record.getSequenceName(),
            start_1based,
            endInclusive_1based);
      }
    }

    byte[] bases = null;

    for (PathPattern pathPattern : cachePatterns) {
      File file = new File(pathPattern.format(md5));
      if (file.exists()) {
        bases = readBytesFromFile(file, start_1based - 1, endInclusive_1based - start_1based + 1);
        return new ReferenceRegion(
            bases, record.getSequenceIndex(), record.getSequenceName(), start_1based);
      }
    }

    { // try to fetch sequence by md5:
      try {
        bases = findBasesRemotelyByMD5(md5);
      } catch (Exception e) {
        if (e instanceof RuntimeException) throw (RuntimeException) e;
        throw new RuntimeException(e);
      }

      if (bases != null) {
        addToMemCache(md5, bases);

        if (!cachePatterns.isEmpty()) addToRefCache(md5, bases);

        return ReferenceRegion.copyRegion(
            bases,
            record.getSequenceIndex(),
            record.getSequenceName(),
            start_1based,
            endInclusive_1based);
      }
    }
    return null;
  }

  private byte[] findBases(SAMSequenceRecord record) {

    String md5 = record.getAttribute(SAMSequenceRecord.MD5_TAG);
    if (null == md5) return null;

    { // check cache by md5:
      byte[] bases = findInMemCache(md5);
      if (bases != null) {
        return bases;
      }
    }

    byte[] bases = null;
    {
      try {
        bases = findInFileCache(md5);
        if (bases != null) {
          addToMemCache(md5, bases);
          return bases;
        }

      } catch (Throwable t) {
        if (t instanceof RuntimeException) throw (RuntimeException) t;
        throw new RuntimeException(t);
      }
    }

    { // try to fetch sequence by md5:
      try {
        bases = findBasesRemotelyByMD5(md5);
        if (bases != null) {
          addToMemCache(md5, bases);

          if (!cachePatterns.isEmpty()) addToRefCache(md5, bases);

          return bases;
        }

      } catch (Throwable t) {
        if (t instanceof RuntimeException) throw (RuntimeException) t;
        throw new RuntimeException(t);
      }
    }

    { // try @SQ:UR file location
      if (record.getAttribute(SAMSequenceRecord.URI_TAG) != null) {
        try {
          ReferenceSequenceFromSeekable s =
              ReferenceSequenceFromSeekable.fromString(
                  record.getAttribute(SAMSequenceRecord.URI_TAG));
          bases = s.getSubsequenceAt(record.getSequenceName(), 1, record.getSequenceLength());
          Utils.upperCase(bases);
          return bases;
        } catch (Throwable t) {
          return null;
        }
      }
    }
    return null;
  }

  /**
   * @param path
   * @return true if the path is a valid URL, false otherwise.
   */
  private static boolean isURL(String path) {
    try {
      URL url = new URL(path);
      return null != url;
    } catch (MalformedURLException e) {
      return false;
    }
  }

  private byte[] loadFromPathWithFileCleanupRetry(File file, String md5) throws IOException {
    RetryTemplate rt = getRetryTemplate(Duration.ofSeconds(1), Duration.ofSeconds(8), 2.0, 3);
    return rt.execute(context -> loadFromPathWithCacheCleanupOnError(file, md5));
  }

  public static RetryTemplate getRetryTemplate(
      Duration minBackoff, Duration maxBackoff, Double multiplier, int attempts) {

    RetryTemplate retryTemplate = new RetryTemplate();
    ExponentialBackOffPolicy backOffPolicy = new ExponentialBackOffPolicy();
    backOffPolicy.setInitialInterval(minBackoff.toMillis());
    backOffPolicy.setMaxInterval(maxBackoff.toMillis());
    backOffPolicy.setMultiplier(multiplier);
    retryTemplate.setBackOffPolicy(backOffPolicy);
    SimpleRetryPolicy retryPolicy = new SimpleRetryPolicy();
    retryPolicy.setMaxAttempts(attempts);
    retryTemplate.setRetryPolicy(retryPolicy);
    return retryTemplate;
  }

  private byte[] loadFromPathWithCacheCleanupOnError(File file, String md5) throws IOException {
    try {
      return loadFromPath(file.getPath(), md5);
    } catch (RuntimeException e) {
      if (e.getMessage().contains("MD5 mismatch for cached file")) {
        log.warn("Deleting corrupt CRAM reference cache file " + file.getAbsolutePath());
        file.delete();
      }
      throw e;
    }
  }

  private byte[] loadFromPath(String path, String md5) throws IOException {
    if (isURL(path)) {
      for (int i = 0; i < downloadTriesBeforeFailing; i++) {
        try {
          URL url = new URL(path);
          URLConnection conn = url.openConnection();
          conn.connect();
          long length = conn.getContentLengthLong();
          InputStream stream = conn.getInputStream();

          if (stream == null || 0 == length) return null;

          try (BufferedInputStream is = new BufferedInputStream(stream)) {
            if (!cachePatterns.isEmpty()) {
              String localPath = addToRefCache(md5, is);
              File file = new File(localPath);
              if (file.length() > Integer.MAX_VALUE) {
                log.warn("The reference sequence is too long: " + url.toExternalForm());
                // throw new RuntimeException( "The reference sequence is too long: " + md5 );
                continue;
              }
              return readBytesFromFile(file, 0, (int) file.length());
            }

            byte[] data =
                length > 0
                    ? InputStreamUtils.readFully(is, (int) length)
                    : InputStreamUtils.readFully(is);

            if (confirmMD5(md5, data)) {
              // sanitise, Internet is a wild place:
              if (Utils.isValidSequence(data, REF_BASES_TO_CHECK_FOR_SANITY)) return data;
              else {
                // reject, it looks like garbage
                log.warn("Downloaded sequence looks suspicous, rejected: " + url.toExternalForm());
                continue;
              }
            }
          }
        } catch (IOException ioe) {
          log.warn(ioe.getClass().getSimpleName() + ". Unable to fetch data from " + path);
          return null;
        }
      }
    } else {
      File file = new File(path);
      if (file.exists()) {
        if (file.length() > Integer.MAX_VALUE)
          throw new RuntimeException("The reference sequence is too long: " + md5);

        byte[] data = readBytesFromFile(file, 0, (int) file.length());

        if (confirmMD5(md5, data)) return data;
        else throw new RuntimeException("MD5 mismatch for cached file: " + file.getAbsolutePath());
      }
    }
    return null;
  }

  private byte[] findBasesRemotelyByMD5(String md5) throws MalformedURLException, IOException {
    long start = System.currentTimeMillis();
    for (PathPattern p : refPatterns) {
      String path = p.format(md5);
      byte[] data = loadFromPath(path, md5);
      if (data == null) continue;

      log.debug(
          String.format(
              "*% 5d Reference found at the location %s sz:%d, total: %d, spent: %d, total spent: %d, attempt %d",
              download_counter.incrementAndGet(),
              path,
              data.length,
              download_sz.addAndGet(data.length),
              System.currentTimeMillis() - start,
              download_spent.addAndGet(System.currentTimeMillis() - start),
              download_map.merge(md5, (Integer) 1, (v1, v2) -> v1 + v2)));
      return data;
    }

    return null;
  }

  private void addToRefCache(String md5, byte[] data) {
    for (PathPattern p : cachePatterns) {
      File cachedFile = new File(p.format(md5));
      if (!cachedFile.exists()) {
        log.debug(String.format("Adding to REF_CACHE: md5=%s, length=%d", md5, data.length));
        cachedFile.getParentFile().mkdirs();
        File tempDir = cachedFile.getParentFile();
        File tmpFile = null;
        try {
          tmpFile = File.createTempFile(md5, ".tmp", tempDir);
        } catch (IOException e) {
          throw new RuntimeException("Error creating temp file in directory : " + tempDir, e);
        }

        try (OutputStream fos = new BufferedOutputStream(new FileOutputStream(tmpFile))) {
          fos.write(data);
          fos.flush();
        } catch (IOException e) {
          throw new RuntimeException(
              "Error creating cached file : " + cachedFile.getAbsolutePath(), e);
        }

        if (!cachedFile.exists()) {
          if (!tmpFile.renameTo(cachedFile)) {
            throw new RuntimeException(
                "'"
                    + tmpFile.getAbsolutePath()
                    + "' rename to '"
                    + cachedFile.getAbsolutePath()
                    + "' failed.");
          }
        } else {
          tmpFile.delete();
        }
      }
    }
  }

  private String addToRefCache(String md5, InputStream stream) {
    for (PathPattern p : cachePatterns) {
      String localPath = p.format(md5);
      File cachedFile = new File(localPath);
      if (!cachedFile.exists()) {
        log.debug(String.format("Adding to REF_CACHE sequence md5=%s", md5));
        cachedFile.getParentFile().mkdirs();
        File tempDir = cachedFile.getParentFile();
        File tmpFile = null;
        try {
          tmpFile = File.createTempFile(md5, ".tmp", tempDir);
        } catch (IOException e) {
          throw new RuntimeException("Error creating temp file in directory : " + tempDir, e);
        }

        try (FileOutputStream fos = new FileOutputStream(tmpFile)) {
          IOUtil.copyStream(stream, fos);
        } catch (IOException e) {
          throw new RuntimeException(
              "Error creating cached file : " + cachedFile.getAbsolutePath(), e);
        }

        if (!cachedFile.exists()) {
          if (!tmpFile.renameTo(cachedFile)) {
            throw new RuntimeException(
                "'"
                    + tmpFile.getAbsolutePath()
                    + "' rename to '"
                    + cachedFile.getAbsolutePath()
                    + "' failed.");
          }
        } else {
          tmpFile.delete();
        }
      }
      return localPath;
    }
    return null;
  }

  private boolean confirmMD5(String md5, byte[] data) {
    String downloadedMD5 = null;
    downloadedMD5 = Utils.calculateMD5String(data);
    if (md5.equals(downloadedMD5)) {
      return true;
    } else {
      String message =
          String.format(
              "Downloaded sequence is corrupt: requested md5=%s, received md5=%s",
              md5, downloadedMD5);
      log.error(message);
      return false;
    }
  }

  public int getDownloadTriesBeforeFailing() {
    return downloadTriesBeforeFailing;
  }

  public void setDownloadTriesBeforeFailing(int downloadTriesBeforeFailing) {
    this.downloadTriesBeforeFailing = downloadTriesBeforeFailing;
  }
}
