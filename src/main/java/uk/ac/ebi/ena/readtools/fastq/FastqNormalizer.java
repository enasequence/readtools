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
package uk.ac.ebi.ena.readtools.fastq;

import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.fastq.AsyncFastqWriter;
import htsjdk.samtools.fastq.BasicFastqWriter;
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.util.FastqQualityFormat;
import java.io.*;
import java.nio.charset.StandardCharsets;
import java.util.*;
import java.util.zip.Deflater;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;
import uk.ac.ebi.ena.readtools.common.reads.CasavaRead;
import uk.ac.ebi.ena.readtools.common.reads.QualityNormalizer;
import uk.ac.ebi.ena.readtools.loader.common.Pair;
import uk.ac.ebi.ena.readtools.loader.common.writer.ReadWriterException;
import uk.ac.ebi.ena.readtools.loader.common.writer.ReadWriterMemoryLimitException;
import uk.ac.ebi.ena.readtools.loader.fastq.PairedFastqWriter;
import uk.ac.ebi.ena.readtools.utils.Utils;

/**
 * FastqNormalizer provides direct FASTQ-to-FASTQ conversion with quality score normalization,
 * uracil base conversion, and optional read name prefixing. This avoids the overhead of BAM
 * intermediate format used in the Fastq2Sam → Sam2Fastq pipeline.
 *
 * <p>Supports both single-end and paired-end FASTQ files. For paired-end, implements full read
 * pairing, sorting, and disk spilling logic matching the behavior of PairedFastqWriter and
 * AbstractPagedReadWriter.
 */
public class FastqNormalizer {

  private static final int OUTPUT_BUFFER_SIZE = 8192;

  /** Result of a paired-end normalization, providing detailed counts. */
  public static class PairedNormalizationResult {
    private final long pairCount;
    private final long orphanCount;
    private final long baseCount;

    public PairedNormalizationResult(long pairCount, long orphanCount, long baseCount) {
      this.pairCount = pairCount;
      this.orphanCount = orphanCount;
      this.baseCount = baseCount;
    }

    /** Number of complete pairs written. */
    public long getPairCount() {
      return pairCount;
    }

    /** Number of orphaned reads (missing mate) written. */
    public long getOrphanCount() {
      return orphanCount;
    }

    /** Total number of individual reads written (pairs * 2 + orphans). */
    public long getTotalReadCount() {
      return pairCount * 2 + orphanCount;
    }

    /** Total number of bases across all written reads. */
    public long getBaseCount() {
      return baseCount;
    }
  }

  /** Result of a single-end normalization, providing detailed counts. */
  public static class SingleNormalizationResult {
    private final long readCount;
    private final long baseCount;

    public SingleNormalizationResult(long readCount, long baseCount) {
      this.readCount = readCount;
      this.baseCount = baseCount;
    }

    /** Number of reads written. */
    public long getReadCount() {
      return readCount;
    }

    /** Total number of bases across all written reads. */
    public long getBaseCount() {
      return baseCount;
    }
  }

  /**
   * Normalizes a single-end FASTQ file.
   *
   * @param inputFastq Path to input FASTQ file (gz/bz2/plain auto-detected)
   * @param outputFastq Path to output FASTQ file (extension determines compression)
   * @param prefix Optional run ID prefix for read names (nullable). When set, read names are
   *     rewritten as "{prefix}.{counter} {originalName}"
   * @param convertUracil If true, converts U bases to T
   * @return Number of reads written
   * @throws IOException If file I/O fails
   */
  public static long normalizeSingleEnd(
      String inputFastq, String outputFastq, String prefix, boolean convertUracil)
      throws IOException {
    return normalizeSingleEndWithStats(inputFastq, outputFastq, prefix, convertUracil)
        .getReadCount();
  }

  /**
   * Normalizes a single-end FASTQ file, returning detailed statistics.
   *
   * @param inputFastq Path to input FASTQ file (gz/bz2/plain auto-detected)
   * @param outputFastq Path to output FASTQ file (extension determines compression)
   * @param prefix Optional run ID prefix for read names (nullable)
   * @param convertUracil If true, converts U bases to T
   * @return Result containing read count and base count
   * @throws IOException If file I/O fails
   */
  public static SingleNormalizationResult normalizeSingleEndWithStats(
      String inputFastq, String outputFastq, String prefix, boolean convertUracil)
      throws IOException {

    // Detect quality format
    FastqQualityFormat format = Utils.detectFastqQualityFormat(inputFastq, null);
    QualityNormalizer normalizer = Utils.getQualityNormalizer(format);

    // Open input and output
    FastqReader reader =
        new FastqReader(
            null,
            new BufferedReader(
                new InputStreamReader(
                    Utils.openFastqInputStream(new File(inputFastq).toPath()),
                    StandardCharsets.UTF_8)));

    AsyncFastqWriter writer =
        new AsyncFastqWriter(
            new BasicFastqWriter(new File(outputFastq)), AsyncFastqWriter.DEFAULT_QUEUE_SIZE);

    long counter = 0;
    long baseCount = 0;

    try {
      for (FastqRecord record : reader) {
        counter++;

        // Process bases
        String bases = record.getReadString();
        if (convertUracil) {
          bases = Utils.replaceUracilBases(bases);
        }
        baseCount += bases.length();

        // Normalize quality scores
        byte[] qualityBytes = record.getBaseQualityString().getBytes(StandardCharsets.UTF_8);
        normalizer.normalize(qualityBytes);
        // phredToFastq converts binary phred scores to ASCII string
        String normalizedQuality = SAMUtils.phredToFastq(qualityBytes);

        // Build read name
        String readName;
        if (prefix != null) {
          readName = prefix + "." + counter + " " + record.getReadName();
        } else {
          readName = record.getReadName();
        }

        // Write record
        writer.write(new FastqRecord(readName, bases, "", normalizedQuality));
      }
    } finally {
      reader.close();
      writer.close();
    }

    return new SingleNormalizationResult(counter, baseCount);
  }

  /**
   * Normalizes paired-end FASTQ files with default thresholds.
   *
   * @param inputFastq1 Path to first mate FASTQ file
   * @param inputFastq2 Path to second mate FASTQ file
   * @param outputFastq1 Path to output first mate FASTQ file
   * @param outputFastq2 Path to output second mate FASTQ file
   * @param prefix Optional run ID prefix (nullable)
   * @param convertUracil If true, converts U bases to T
   * @param tempDir Directory for temporary spill files
   * @return Result containing pair count, orphan count, and total read count
   * @throws IOException If file I/O fails
   */
  public static PairedNormalizationResult normalizePairedEnd(
      String inputFastq1,
      String inputFastq2,
      String outputFastq1,
      String outputFastq2,
      String prefix,
      boolean convertUracil,
      File tempDir)
      throws IOException {
    return normalizePairedEnd(
        inputFastq1,
        inputFastq2,
        outputFastq1,
        outputFastq2,
        prefix,
        convertUracil,
        tempDir,
        100_000, // 100K reads in memory
        4L * 1024L * 1024L * 1024L, // 4 GB memory threshold
        10L * 1024L * 1024L * 1024L); // 10 GB disk abandon limit
  }

  /**
   * Normalizes paired-end FASTQ files with full pairing/sorting/spilling support.
   *
   * @param inputFastq1 Path to first mate FASTQ file
   * @param inputFastq2 Path to second mate FASTQ file
   * @param outputFastq1 Path to output first mate FASTQ file
   * @param outputFastq2 Path to output second mate FASTQ file
   * @param prefix Optional run ID prefix (nullable)
   * @param convertUracil If true, converts U bases to T
   * @param tempDir Directory for temporary spill files
   * @param spillPageSize Maximum number of reads to keep in memory before spilling
   * @param spillPageSizeBytes Maximum memory usage in bytes before spilling
   * @param spillAbandonLimitBytes Maximum total spilled bytes before aborting
   * @return Result containing pair count, orphan count, and total read count
   * @throws IOException If file I/O fails
   */
  public static PairedNormalizationResult normalizePairedEnd(
      String inputFastq1,
      String inputFastq2,
      String outputFastq1,
      String outputFastq2,
      String prefix,
      boolean convertUracil,
      File tempDir,
      int spillPageSize,
      long spillPageSizeBytes,
      long spillAbandonLimitBytes)
      throws IOException {

    PairedNormalizer normalizer =
        new PairedNormalizer(
            inputFastq1,
            inputFastq2,
            outputFastq1,
            outputFastq2,
            prefix,
            convertUracil,
            tempDir,
            spillPageSize,
            spillPageSizeBytes,
            spillAbandonLimitBytes);

    return normalizer.normalize();
  }

  /** Helper class for paired-end normalization with buffering and spilling. */
  private static class PairedNormalizer {
    private final String inputFastq1;
    private final String inputFastq2;
    private final String outputFastq1;
    private final String outputFastq2;
    private final String prefix;
    private final boolean convertUracil;
    private final File tempDir;
    private final int spillPageSize;
    private final long spillPageSizeBytes;
    private final long spillAbandonLimitBytes;

    private QualityNormalizer qualityNormalizer;
    private Map<String, List<NormalizedRead>> pairMap;
    private long totalBytesInMemory;
    private long totalSpilledBytes;
    private List<File> spillFiles;
    private Integer index1;
    private Integer index2;

    public PairedNormalizer(
        String inputFastq1,
        String inputFastq2,
        String outputFastq1,
        String outputFastq2,
        String prefix,
        boolean convertUracil,
        File tempDir,
        int spillPageSize,
        long spillPageSizeBytes,
        long spillAbandonLimitBytes) {
      this.inputFastq1 = inputFastq1;
      this.inputFastq2 = inputFastq2;
      this.outputFastq1 = outputFastq1;
      this.outputFastq2 = outputFastq2;
      this.prefix = prefix;
      this.convertUracil = convertUracil;
      this.tempDir = tempDir;
      this.spillPageSize = spillPageSize;
      this.spillPageSizeBytes = spillPageSizeBytes;
      this.spillAbandonLimitBytes = spillAbandonLimitBytes;
      this.pairMap = new HashMap<>(spillPageSize);
      this.totalBytesInMemory = 0;
      this.totalSpilledBytes = 0;
      this.spillFiles = new ArrayList<>();
      this.index1 = null;
      this.index2 = null;
    }

    public PairedNormalizationResult normalize() throws IOException {
      // Detect quality format
      FastqQualityFormat format = Utils.detectFastqQualityFormat(inputFastq1, inputFastq2);
      qualityNormalizer = Utils.getQualityNormalizer(format);

      // Process both input files
      processInputFiles();

      // Open output writers shared across in-memory write and spill processing
      AsyncFastqWriter writer1 =
          new AsyncFastqWriter(
              new BasicFastqWriter(new File(outputFastq1)), AsyncFastqWriter.DEFAULT_QUEUE_SIZE);
      AsyncFastqWriter writer2 =
          new AsyncFastqWriter(
              new BasicFastqWriter(new File(outputFastq2)), AsyncFastqWriter.DEFAULT_QUEUE_SIZE);

      long counter = 0;
      long pairCount = 0;
      long orphanCount = 0;
      long baseCount = 0;

      try {
        if (spillFiles.isEmpty()) {
          // No spilling occurred — write everything from memory
          WriteCounts counts = writeFromMemory(writer1, writer2, counter);
          counter = counts.counter;
          pairCount += counts.pairCount;
          orphanCount += counts.orphanCount;
          baseCount += counts.baseCount;
        } else {
          // Spilling occurred — the residual in-memory data is the base for reassembly.
          // processSpillFiles() will use pairMap as the starting point, matching the
          // pattern in AbstractPagedReadWriter.cascadeErrors().
          WriteCounts spillCounts = processSpillFiles(writer1, writer2, counter);
          counter = spillCounts.counter;
          pairCount += spillCounts.pairCount;
          orphanCount += spillCounts.orphanCount;
          baseCount += spillCounts.baseCount;
        }
      } finally {
        writer1.close();
        writer2.close();

        // Clean up spill files
        for (File spillFile : spillFiles) {
          spillFile.delete();
        }
      }

      return new PairedNormalizationResult(pairCount, orphanCount, baseCount);
    }

    private void processInputFiles() throws IOException {
      FastqReader reader1 =
          new FastqReader(
              null,
              new BufferedReader(
                  new InputStreamReader(
                      Utils.openFastqInputStream(new File(inputFastq1).toPath()),
                      StandardCharsets.UTF_8)));

      FastqReader reader2 =
          new FastqReader(
              null,
              new BufferedReader(
                  new InputStreamReader(
                      Utils.openFastqInputStream(new File(inputFastq2).toPath()),
                      StandardCharsets.UTF_8)));

      try {
        Iterator<FastqRecord> iter1 = reader1.iterator();
        Iterator<FastqRecord> iter2 = reader2.iterator();

        // Read both files round-robin, tolerating different lengths.
        // This matches MultiFastqConverter behavior: when one file ends,
        // continue reading the other.
        while (iter1.hasNext() || iter2.hasNext()) {
          if (iter1.hasNext()) {
            processRecord(iter1.next());
          }
          if (iter2.hasNext()) {
            processRecord(iter2.next());
          }
        }
      } finally {
        reader1.close();
        reader2.close();
      }
    }

    private void processRecord(FastqRecord record) throws IOException {
      try {
        // Normalize bases
        String bases = record.getReadString();
        if (convertUracil) {
          bases = Utils.replaceUracilBases(bases);
        }

        // Normalize quality
        byte[] qualityBytes = record.getBaseQualityString().getBytes(StandardCharsets.UTF_8);
        qualityNormalizer.normalize(qualityBytes);
        String normalizedQuality = SAMUtils.phredToFastq(qualityBytes);

        // Extract read key and pair number
        String readKey = PairedFastqWriter.getReadKey(record.getReadName());
        String pairNumberStr = PairedFastqWriter.getPairNumber(record.getReadName());
        int pairNumber = Integer.parseInt(pairNumberStr);

        // Track pair indices
        if (index1 == null) {
          index1 = pairNumber;
        } else if (index2 == null && pairNumber != index1) {
          index2 = pairNumber;
        }

        // Validate pair number
        if (pairNumber != index1 && pairNumber != index2) {
          throw new IOException(
              "Unexpected read pair number: "
                  + pairNumber
                  + "; pair numbers "
                  + index1
                  + " and "
                  + index2
                  + " were found previously");
        }

        // Create normalized read
        NormalizedRead normRead =
            new NormalizedRead(record.getReadName(), bases, normalizedQuality, pairNumber);

        // Add to buffer
        List<NormalizedRead> readList = pairMap.get(readKey);
        if (readList == null) {
          readList = new ArrayList<>(2);
          readList.add(null);
          readList.add(null);
          pairMap.put(readKey, readList);
        }

        int mappedIndex = (pairNumber == index1) ? 0 : 1;
        if (readList.get(mappedIndex) != null) {
          throw new IOException("Got same spot twice: " + record.getReadName());
        }
        readList.set(mappedIndex, normRead);

        // Update memory tracking
        totalBytesInMemory +=
            normRead.readName.length()
                + normRead.bases.length()
                + normRead.qualities.length()
                + 100;

        // Check if we need to spill
        if (pairMap.size() >= spillPageSize || totalBytesInMemory >= spillPageSizeBytes) {
          spillToDisk();
        }

      } catch (ReadWriterException e) {
        throw new IOException("Failed to extract read key/pair number: " + e.getMessage(), e);
      }
    }

    private void spillToDisk() throws IOException {
      if (spillAbandonLimitBytes > 0 && totalSpilledBytes >= spillAbandonLimitBytes) {
        throw new ReadWriterMemoryLimitException(
            "Temp memory limit " + spillAbandonLimitBytes + " bytes reached");
      }

      File spillFile = createTempFile();
      spillFiles.add(spillFile);

      ObjectOutputStream oos = openOutputStream(spillFile);
      try {
        for (Map.Entry<String, List<NormalizedRead>> entry : pairMap.entrySet()) {
          oos.writeObject(new Pair<>(entry.getKey(), entry.getValue()));
          oos.reset();
        }
      } finally {
        oos.close();
      }

      totalSpilledBytes += totalBytesInMemory;
      pairMap.clear();
      totalBytesInMemory = 0;
    }

    /** Writes pairs and orphans from the in-memory pairMap, then clears it. */
    private WriteCounts writeFromMemory(
        AsyncFastqWriter writer1, AsyncFastqWriter writer2, long counter) {
      // Sort keys lexicographically to match BAM queryname sort order
      List<String> sortedKeys = new ArrayList<>(pairMap.keySet());
      Collections.sort(sortedKeys);

      long pairCount = 0;
      long orphanCount = 0;
      long baseCount = 0;

      for (String key : sortedKeys) {
        List<NormalizedRead> reads = pairMap.get(key);
        counter++;

        if (reads.get(0) != null && reads.get(1) != null) {
          writePair(writer1, writer2, reads, counter);
          baseCount += reads.get(0).bases.length() + reads.get(1).bases.length();
          pairCount++;
        } else {
          NormalizedRead orphan = reads.get(0) != null ? reads.get(0) : reads.get(1);
          writeOrphan(writer1, orphan, counter);
          baseCount += orphan.bases.length();
          orphanCount++;
        }
      }

      pairMap.clear();
      totalBytesInMemory = 0;

      return new WriteCounts(counter, pairCount, orphanCount, baseCount);
    }

    /**
     * Multi-generation spill file reassembly, following the pattern from
     * AbstractPagedReadWriter.cascadeErrors().
     *
     * <p>Algorithm: load one spill file into pairMap as the base. Stream each remaining spill file
     * entry-by-entry: if the key exists in pairMap, merge the reads (completing a pair or adding to
     * the slot); if not, re-spill to a new generation file. After all files in one generation are
     * processed, write completed pairs and orphans from pairMap, then repeat with any new
     * generation files.
     *
     * @return WriteCounts with accumulated counter, pairCount, orphanCount, baseCount
     */
    private WriteCounts processSpillFiles(
        AsyncFastqWriter writer1, AsyncFastqWriter writer2, long counter) throws IOException {
      long pairCount = 0;
      long orphanCount = 0;
      long baseCount = 0;

      int i = 0;
      do {
        // If pairMap is empty, load the first unprocessed spill file as base
        if (pairMap.isEmpty()) {
          loadSpillFile(spillFiles.get(i++));
        }

        // Snapshot current file count — files added during this loop are next generation
        int generation = spillFiles.size();

        // Stream remaining files in this generation against the in-memory base
        for (int j = i; j < generation; j++) {
          ObjectInputStream ois = null;
          ObjectOutputStream oos = null;

          try {
            ois = openInputStream(spillFiles.get(j));

            for (; ; ) {
              @SuppressWarnings("unchecked")
              Pair<String, List<NormalizedRead>> entry =
                  (Pair<String, List<NormalizedRead>>) ois.readObject();

              if (pairMap.containsKey(entry.key)) {
                // Merge: fill in the missing slot(s)
                List<NormalizedRead> existing = pairMap.get(entry.key);
                for (int k = 0; k < entry.value.size(); k++) {
                  if (entry.value.get(k) != null && existing.get(k) == null) {
                    existing.set(k, entry.value.get(k));
                  }
                }
              } else {
                // No match in memory — re-spill to next generation
                if (oos == null) {
                  File nextGenFile = createTempFile();
                  spillFiles.add(nextGenFile);
                  oos = openOutputStream(nextGenFile);
                }
                oos.writeObject(entry);
                entry.value = null;
                oos.reset();
              }
            }
          } catch (EOFException eof) {
            // Normal end of spill file
            if (ois != null) {
              ois.close();
            }
          } catch (ClassNotFoundException e) {
            throw new IOException("Failed to deserialize spill file entry", e);
          }

          if (oos != null) {
            oos.close();
          }
        }

        // Write completed pairs and remaining orphans from this generation
        WriteCounts counts = writeFromMemory(writer1, writer2, counter);
        counter = counts.counter;
        pairCount += counts.pairCount;
        orphanCount += counts.orphanCount;
        baseCount += counts.baseCount;

        // Advance to next generation
        i = generation;
      } while (i < spillFiles.size());

      return new WriteCounts(counter, pairCount, orphanCount, baseCount);
    }

    /** Loads a spill file into pairMap. */
    @SuppressWarnings("unchecked")
    private void loadSpillFile(File file) throws IOException {
      pairMap.clear();
      totalBytesInMemory = 0;

      try (ObjectInputStream ois = openInputStream(file)) {
        for (; ; ) {
          Pair<String, List<NormalizedRead>> entry =
              (Pair<String, List<NormalizedRead>>) ois.readObject();
          pairMap.put(entry.key, entry.value);
          for (NormalizedRead read : entry.value) {
            if (read != null) {
              totalBytesInMemory +=
                  read.readName.length() + read.bases.length() + read.qualities.length() + 100;
            }
          }
        }
      } catch (EOFException eof) {
        // Normal end of file
      } catch (ClassNotFoundException e) {
        throw new IOException("Failed to deserialize spill file", e);
      }
    }

    private void writePair(
        AsyncFastqWriter writer1,
        AsyncFastqWriter writer2,
        List<NormalizedRead> reads,
        long counter) {

      // Sort by pair number
      reads.sort(Comparator.comparingInt(r -> r.pairNumber));

      NormalizedRead read1 = reads.get(0);
      NormalizedRead read2 = reads.get(1);

      // Casava 1.8 read names already contain the pair number in the metadata
      // (e.g. "inst:run:fc:lane:tile:x:y 1:N:0:barcode"), so don't append /1 /2.
      boolean casava = CasavaRead.getBaseNameOrNull(read1.readName) != null;

      String name1, name2;
      if (casava) {
        if (prefix != null) {
          name1 = prefix + "." + counter + " " + read1.readName;
          name2 = prefix + "." + counter + " " + read2.readName;
        } else {
          name1 = read1.readName;
          name2 = read2.readName;
        }
      } else {
        String baseName1 = stripPairSuffix(read1.readName);
        String baseName2 = stripPairSuffix(read2.readName);
        if (prefix != null) {
          name1 = prefix + "." + counter + " " + baseName1 + "/1";
          name2 = prefix + "." + counter + " " + baseName2 + "/2";
        } else {
          name1 = baseName1 + "/1";
          name2 = baseName2 + "/2";
        }
      }

      writer1.write(new FastqRecord(name1, read1.bases, "", read1.qualities));
      writer2.write(new FastqRecord(name2, read2.bases, "", read2.qualities));
    }

    private void writeOrphan(AsyncFastqWriter writer, NormalizedRead read, long counter) {
      // Casava read names are kept as-is; non-Casava get pair suffix stripped.
      boolean casava = CasavaRead.getBaseNameOrNull(read.readName) != null;
      String baseName = casava ? read.readName : stripPairSuffix(read.readName);

      String name;
      if (prefix != null) {
        name = prefix + "." + counter + " " + baseName;
      } else {
        name = baseName;
      }

      writer.write(new FastqRecord(name, read.bases, "", read.qualities));
    }

    private static String stripPairSuffix(String readName) {
      try {
        return PairedFastqWriter.getReadKey(readName);
      } catch (ReadWriterException e) {
        // If the name doesn't match any known pattern, return as-is.
        return readName;
      }
    }

    private File createTempFile() throws IOException {
      String prefix =
          String.format(
              "FASTQ_NORM_THREAD_%d_TIME_%d_FILE_",
              Thread.currentThread().getId(), System.currentTimeMillis());
      String suffix = String.format("_PAGE_%d", spillFiles.size());
      File tmpFile = File.createTempFile(prefix, suffix, tempDir);
      tmpFile.deleteOnExit();
      return tmpFile;
    }

    private ObjectOutputStream openOutputStream(File file) throws IOException {
      return new ObjectOutputStream(
          new BufferedOutputStream(
              new GZIPOutputStream(new BufferedOutputStream(new FileOutputStream(file))) {
                {
                  def.setLevel(Deflater.BEST_SPEED);
                }
              },
              OUTPUT_BUFFER_SIZE));
    }

    private ObjectInputStream openInputStream(File file) throws IOException {
      return new ObjectInputStream(
          new BufferedInputStream(
              new GZIPInputStream(new BufferedInputStream(new FileInputStream(file)))));
    }
  }

  /** Accumulated counts from a batch of writes. */
  private static class WriteCounts {
    final long counter;
    final long pairCount;
    final long orphanCount;
    final long baseCount;

    WriteCounts(long counter, long pairCount, long orphanCount, long baseCount) {
      this.counter = counter;
      this.pairCount = pairCount;
      this.orphanCount = orphanCount;
      this.baseCount = baseCount;
    }
  }

  /** Represents a normalized FASTQ read ready for output. */
  private static class NormalizedRead implements Serializable {
    private static final long serialVersionUID = 1L;

    final String readName;
    final String bases; // After uracil conversion
    final String qualities; // After normalization to Phred+33
    final int pairNumber;

    public NormalizedRead(String readName, String bases, String qualities, int pairNumber) {
      this.readName = readName;
      this.bases = bases;
      this.qualities = qualities;
      this.pairNumber = pairNumber;
    }
  }
}
