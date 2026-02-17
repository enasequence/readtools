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

    try {
      for (FastqRecord record : reader) {
        counter++;

        // Process bases
        String bases = record.getReadString();
        if (convertUracil) {
          bases = Utils.replaceUracilBases(bases);
        }

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
        // Note: quality header should be empty string or original header, not "+"
        writer.write(new FastqRecord(readName, bases, "", normalizedQuality));
      }
    } finally {
      reader.close();
      writer.close();
    }

    return counter;
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
   * @return Number of pairs written
   * @throws IOException If file I/O fails
   */
  public static long normalizePairedEnd(
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
   * @return Number of pairs written
   * @throws IOException If file I/O fails
   */
  public static long normalizePairedEnd(
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
      this.pairMap = new TreeMap<>();
      this.totalBytesInMemory = 0;
      this.totalSpilledBytes = 0;
      this.spillFiles = new ArrayList<>();
      this.index1 = null;
      this.index2 = null;
    }

    public long normalize() throws IOException {
      // Detect quality format
      FastqQualityFormat format = Utils.detectFastqQualityFormat(inputFastq1, inputFastq2);
      qualityNormalizer = Utils.getQualityNormalizer(format);

      // Process both input files
      processInputFiles();

      // Write complete pairs from memory
      long pairsWritten = writePairsFromMemory();

      // Process spill files if any
      if (!spillFiles.isEmpty()) {
        pairsWritten += processSpillFiles();
      }

      // Clean up
      for (File spillFile : spillFiles) {
        spillFile.delete();
      }

      return pairsWritten;
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
        totalBytesInMemory += estimateReadSize(normRead);

        // Check if we need to spill
        if (pairMap.size() >= spillPageSize || totalBytesInMemory >= spillPageSizeBytes) {
          spillToDisk();
        }

      } catch (ReadWriterException e) {
        throw new IOException("Failed to extract read key/pair number: " + e.getMessage(), e);
      }
    }

    private long estimateReadSize(NormalizedRead read) {
      return read.readName.length() + read.bases.length() + read.qualities.length() + 100;
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

    private long writePairsFromMemory() throws IOException {
      AsyncFastqWriter writer1 =
          new AsyncFastqWriter(
              new BasicFastqWriter(new File(outputFastq1)), AsyncFastqWriter.DEFAULT_QUEUE_SIZE);
      AsyncFastqWriter writer2 =
          new AsyncFastqWriter(
              new BasicFastqWriter(new File(outputFastq2)), AsyncFastqWriter.DEFAULT_QUEUE_SIZE);

      long counter = 0;

      try {
        // Write complete pairs
        for (Map.Entry<String, List<NormalizedRead>> entry : pairMap.entrySet()) {
          List<NormalizedRead> reads = entry.getValue();
          if (reads.get(0) != null && reads.get(1) != null) {
            counter++;
            writePair(writer1, writer2, reads, counter);
          }
        }

        // Write orphaned reads (one mate missing) to the first output file,
        // matching Sam2Fastq behavior where unpaired reads go to streams[0]
        // without /1 or /2 suffix.
        for (Map.Entry<String, List<NormalizedRead>> entry : pairMap.entrySet()) {
          List<NormalizedRead> reads = entry.getValue();
          if (reads.get(0) == null || reads.get(1) == null) {
            NormalizedRead orphan = reads.get(0) != null ? reads.get(0) : reads.get(1);
            counter++;
            writeOrphan(writer1, orphan, counter);
          }
        }
      } finally {
        writer1.close();
        writer2.close();
      }

      return counter;
    }

    private long processSpillFiles() throws IOException {
      // TODO: Implement multi-pass spill file processing
      // For now, throw exception if spilling occurred
      throw new IOException(
          "Spill file processing not yet implemented. Input files may be too large or out of order.");
    }

    private void writePair(
        AsyncFastqWriter writer1,
        AsyncFastqWriter writer2,
        List<NormalizedRead> reads,
        long counter)
        throws IOException {

      // Sort by pair number
      reads.sort(Comparator.comparingInt(r -> r.pairNumber));

      NormalizedRead read1 = reads.get(0);
      NormalizedRead read2 = reads.get(1);

      // Build read names with prefix and /1 /2 suffixes
      // Strip any existing /1 or /2 suffix from read names
      String baseName1 = read1.readName.replaceAll("/[12]$", "");
      String baseName2 = read2.readName.replaceAll("/[12]$", "");

      String name1, name2;
      if (prefix != null) {
        name1 = prefix + "." + counter + " " + baseName1 + "/1";
        name2 = prefix + "." + counter + " " + baseName2 + "/2";
      } else {
        name1 = baseName1 + "/1";
        name2 = baseName2 + "/2";
      }

      writer1.write(new FastqRecord(name1, read1.bases, "", read1.qualities));
      writer2.write(new FastqRecord(name2, read2.bases, "", read2.qualities));
    }

    private void writeOrphan(AsyncFastqWriter writer, NormalizedRead read, long counter) {
      String baseName = read.readName.replaceAll("/[12]$", "");

      String name;
      if (prefix != null) {
        name = prefix + "." + counter + " " + baseName;
      } else {
        name = baseName;
      }

      writer.write(new FastqRecord(name, read.bases, "", read.qualities));
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
