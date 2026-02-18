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

import static org.junit.Assert.*;

import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;
import uk.ac.ebi.ena.readtools.fastq.ena.Fastq2Sam;
import uk.ac.ebi.ena.readtools.loader.common.FileCompression;
import uk.ac.ebi.ena.readtools.sam.Sam2Fastq;
import uk.ac.ebi.ena.readtools.utils.Utils;

/**
 * Compares the existing Fastq2Sam → Sam2Fastq pipeline output against FastqNormalizer output for
 * the same input. Both paths should produce equivalent results for bases, qualities, read counts,
 * and pairing.
 *
 * <p>Known intentional difference: Casava 1.8 read names. The pipeline strips the Casava metadata
 * tail (barcode, filter, control) through the BAM round-trip, while FastqNormalizer preserves it.
 */
public class RoundTripComparisonTest {

  @Rule public TemporaryFolder tempFolder = new TemporaryFolder();

  // ---- Non-Casava paired tests ----

  /** Standard /1 /2 paired reads: both paths should produce identical output. */
  @Test
  public void testSlashSeparatorPaired() throws Exception {
    String[][] reads = {
      {"@READAAA/1", "ACGTACGT", "IIIIIIII", "@READAAA/2", "TTAATTAA", "KKKKKKKK"},
      {"@READBBB/1", "GGCCGGCC", "JJJJJJJJ", "@READBBB/2", "CCGGCCGG", "LLLLLLLL"},
      {"@READCCC/1", "AACCTTGG", "HHHHHHHH", "@READCCC/2", "GGTTCCAA", "MMMMMMMM"},
    };

    RoundTripResult result = runBothPaths(reads, "ERR111", false);

    assertEquals("Read count mismatch", result.pipelinePairs.size(), result.normalizerPairs.size());

    for (String key : result.pipelinePairs.keySet()) {
      assertTrue(
          "Missing key in normalizer output: " + key, result.normalizerPairs.containsKey(key));
      ReadPair pipeline = result.pipelinePairs.get(key);
      ReadPair normalizer = result.normalizerPairs.get(key);

      assertEquals("Bases mismatch for " + key + " /1", pipeline.bases1, normalizer.bases1);
      assertEquals("Bases mismatch for " + key + " /2", pipeline.bases2, normalizer.bases2);
      assertEquals("Quality mismatch for " + key + " /1", pipeline.qual1, normalizer.qual1);
      assertEquals("Quality mismatch for " + key + " /2", pipeline.qual2, normalizer.qual2);
    }
  }

  /** Dot-separated pair index: READ.1 / READ.2 */
  @Test
  public void testDotSeparatorPaired() throws Exception {
    String[][] reads = {
      {"@READX.1", "ACGTACGT", "IIIIIIII", "@READX.2", "TTAATTAA", "KKKKKKKK"},
      {"@READY.1", "GGCCGGCC", "JJJJJJJJ", "@READY.2", "CCGGCCGG", "LLLLLLLL"},
    };

    RoundTripResult result = runBothPaths(reads, "ERR222", false);

    assertEquals(result.pipelinePairs.size(), result.normalizerPairs.size());

    for (String key : result.pipelinePairs.keySet()) {
      assertTrue("Missing key in normalizer: " + key, result.normalizerPairs.containsKey(key));
      ReadPair pipeline = result.pipelinePairs.get(key);
      ReadPair normalizer = result.normalizerPairs.get(key);

      assertEquals(pipeline.bases1, normalizer.bases1);
      assertEquals(pipeline.bases2, normalizer.bases2);
      assertEquals(pipeline.qual1, normalizer.qual1);
      assertEquals(pipeline.qual2, normalizer.qual2);
    }
  }

  /** Underscore-separated pair index: READ_1 / READ_2 */
  @Test
  public void testUnderscoreSeparatorPaired() throws Exception {
    String[][] reads = {
      {"@MYREAD_1", "ACGTACGT", "IIIIIIII", "@MYREAD_2", "TTAATTAA", "KKKKKKKK"},
    };

    RoundTripResult result = runBothPaths(reads, "ERR333", false);

    assertEquals(result.pipelinePairs.size(), result.normalizerPairs.size());

    for (String key : result.pipelinePairs.keySet()) {
      assertTrue(result.normalizerPairs.containsKey(key));
      ReadPair pipeline = result.pipelinePairs.get(key);
      ReadPair normalizer = result.normalizerPairs.get(key);

      assertEquals(pipeline.bases1, normalizer.bases1);
      assertEquals(pipeline.bases2, normalizer.bases2);
    }
  }

  /** Uracil conversion: both paths should produce identical T-converted bases. */
  @Test
  public void testUracilConversion() throws Exception {
    String[][] reads = {
      {"@READ1/1", "AUGCUACG", "IIIIIIII", "@READ1/2", "UUAAuuaa", "KKKKKKKK"},
    };

    RoundTripResult result = runBothPaths(reads, "ERR444", true);

    assertEquals(result.pipelinePairs.size(), result.normalizerPairs.size());

    for (String key : result.pipelinePairs.keySet()) {
      ReadPair pipeline = result.pipelinePairs.get(key);
      ReadPair normalizer = result.normalizerPairs.get(key);

      // Both should have U→T conversion applied
      assertFalse("Pipeline bases still contain U", pipeline.bases1.toUpperCase().contains("U"));
      assertFalse(
          "Normalizer bases still contain U", normalizer.bases1.toUpperCase().contains("U"));
      // Compare case-insensitively: pipeline uppercases through BAM,
      // normalizer preserves original casing of non-uracil bases.
      assertEquals(pipeline.bases1.toUpperCase(), normalizer.bases1.toUpperCase());
      assertEquals(pipeline.bases2.toUpperCase(), normalizer.bases2.toUpperCase());
    }
  }

  // ---- Casava 1.8 paired test ----

  /**
   * Casava 1.8 format: bases and qualities should match. Read names will differ — the pipeline
   * strips the Casava tail through BAM, while FastqNormalizer preserves it.
   */
  @Test
  public void testCasavaFormat() throws Exception {
    Path inputFile1 = tempFolder.newFile("casava_1.fastq").toPath();
    Path inputFile2 = tempFolder.newFile("casava_2.fastq").toPath();

    Files.write(
        inputFile1,
        ("@A00500:310:HG3JFDRXY:1:2101:1000:2000 1:N:0:ATCACG\n"
                + "ACGTACGT\n+\nIIIIIIII\n"
                + "@A00500:310:HG3JFDRXY:1:2101:1000:3000 1:N:0:TTAGGC\n"
                + "GGCCGGCC\n+\nJJJJJJJJ\n")
            .getBytes(StandardCharsets.UTF_8));

    Files.write(
        inputFile2,
        ("@A00500:310:HG3JFDRXY:1:2101:1000:2000 2:N:0:ATCACG\n"
                + "TTAATTAA\n+\nKKKKKKKK\n"
                + "@A00500:310:HG3JFDRXY:1:2101:1000:3000 2:N:0:TTAGGC\n"
                + "CCGGCCGG\n+\nLLLLLLLL\n")
            .getBytes(StandardCharsets.UTF_8));

    String prefix = "ERR555";

    // Run pipeline
    List<FastqRecord> pipelineOut1 = new ArrayList<>();
    List<FastqRecord> pipelineOut2 = new ArrayList<>();
    runPipeline(inputFile1, inputFile2, prefix, false, pipelineOut1, pipelineOut2);

    // Run normalizer
    Path normOut1 = tempFolder.newFile("norm_1.fastq").toPath();
    Path normOut2 = tempFolder.newFile("norm_2.fastq").toPath();
    FastqNormalizer.normalizePairedEnd(
        inputFile1.toString(),
        inputFile2.toString(),
        normOut1.toString(),
        normOut2.toString(),
        prefix,
        false,
        tempFolder.getRoot());

    List<FastqRecord> normRecords1 = readFastq(normOut1);
    List<FastqRecord> normRecords2 = readFastq(normOut2);

    assertEquals("Pair count mismatch", pipelineOut1.size(), normRecords1.size());
    assertEquals("Pair count mismatch", pipelineOut2.size(), normRecords2.size());

    // Sort both by bases (since ordering should match but let's be safe)
    pipelineOut1.sort(Comparator.comparing(FastqRecord::getReadString));
    pipelineOut2.sort(Comparator.comparing(FastqRecord::getReadString));
    normRecords1.sort(Comparator.comparing(FastqRecord::getReadString));
    normRecords2.sort(Comparator.comparing(FastqRecord::getReadString));

    for (int i = 0; i < pipelineOut1.size(); i++) {
      // Bases and qualities must match
      assertEquals(
          "Bases mismatch at index " + i,
          pipelineOut1.get(i).getReadString(),
          normRecords1.get(i).getReadString());
      assertEquals(
          "Quality mismatch at index " + i,
          pipelineOut1.get(i).getBaseQualityString(),
          normRecords1.get(i).getBaseQualityString());
      assertEquals(pipelineOut2.get(i).getReadString(), normRecords2.get(i).getReadString());
      assertEquals(
          pipelineOut2.get(i).getBaseQualityString(), normRecords2.get(i).getBaseQualityString());

      // Read names differ intentionally: pipeline has just instrument name,
      // normalizer preserves the full Casava tail.
      String pipelineName = pipelineOut1.get(i).getReadName();
      String normName = normRecords1.get(i).getReadName();
      assertNotEquals(
          "Casava read names should differ (pipeline strips metadata)", pipelineName, normName);

      // But both should start with the same prefix.counter
      String pipelinePrefix = pipelineName.split(" ")[0];
      String normPrefix = normName.split(" ")[0];
      assertEquals("Prefix.counter should match", pipelinePrefix, normPrefix);
    }
  }

  // ---- Orphan handling test ----

  /**
   * With orphans (no mate), both paths should produce the same bases/qualities for orphaned reads.
   * Note: the pipeline stops reading when the shorter file hits EOF, while FastqNormalizer
   * continues reading, so we use equal-length files with deliberate non-matching keys.
   */
  @Test
  public void testWithOrphans() throws Exception {
    Path inputFile1 = tempFolder.newFile("orphan_1.fastq").toPath();
    Path inputFile2 = tempFolder.newFile("orphan_2.fastq").toPath();

    // Both files have 3 reads. A and C pair, B in file1 and X in file2 are orphans.
    Files.write(
        inputFile1,
        ("@READA/1\nACGTACGT\n+\nIIIIIIII\n"
                + "@READB/1\nGGCCGGCC\n+\nJJJJJJJJ\n"
                + "@READC/1\nAACCTTGG\n+\nHHHHHHHH\n")
            .getBytes(StandardCharsets.UTF_8));

    Files.write(
        inputFile2,
        ("@READA/2\nTTAATTAA\n+\nKKKKKKKK\n"
                + "@XORPHAN/2\nCCGGCCGG\n+\nLLLLLLLL\n"
                + "@READC/2\nGGTTCCAA\n+\nMMMMMMMM\n")
            .getBytes(StandardCharsets.UTF_8));

    String prefix = "ERR666";

    // Run pipeline
    List<FastqRecord> pipelineOut1 = new ArrayList<>();
    List<FastqRecord> pipelineOut2 = new ArrayList<>();
    runPipeline(inputFile1, inputFile2, prefix, false, pipelineOut1, pipelineOut2);

    // Run normalizer
    Path normOut1 = tempFolder.newFile("norm_1.fastq").toPath();
    Path normOut2 = tempFolder.newFile("norm_2.fastq").toPath();
    FastqNormalizer.normalizePairedEnd(
        inputFile1.toString(),
        inputFile2.toString(),
        normOut1.toString(),
        normOut2.toString(),
        prefix,
        false,
        tempFolder.getRoot());

    List<FastqRecord> normRecords1 = readFastq(normOut1);
    List<FastqRecord> normRecords2 = readFastq(normOut2);

    // Both should have the same total read count
    int pipelineTotal = pipelineOut1.size() + pipelineOut2.size();
    int normTotal = normRecords1.size() + normRecords2.size();
    assertEquals("Total read count mismatch", pipelineTotal, normTotal);

    // Collect all bases from both paths and compare as sets
    Set<String> pipelineBases = new HashSet<>();
    for (FastqRecord r : pipelineOut1) pipelineBases.add(r.getReadString());
    for (FastqRecord r : pipelineOut2) pipelineBases.add(r.getReadString());

    Set<String> normBases = new HashSet<>();
    for (FastqRecord r : normRecords1) normBases.add(r.getReadString());
    for (FastqRecord r : normRecords2) normBases.add(r.getReadString());

    assertEquals("Base content should match across both paths", pipelineBases, normBases);
  }

  // ---- Real resource file test ----

  /** Round-trip comparison using actual test resource files. */
  @Test
  public void testRealResourceFiles() throws Exception {
    Path inputFile1 = new File("src/test/resources/2fastq/28239_1822_1.fastq").toPath();
    Path inputFile2 = new File("src/test/resources/2fastq/28239_1822_2.fastq").toPath();

    String prefix = "ERR777";

    // Run pipeline
    List<FastqRecord> pipelineOut1 = new ArrayList<>();
    List<FastqRecord> pipelineOut2 = new ArrayList<>();
    runPipeline(inputFile1, inputFile2, prefix, false, pipelineOut1, pipelineOut2);

    // Run normalizer
    Path normOut1 = tempFolder.newFile("norm_1.fastq").toPath();
    Path normOut2 = tempFolder.newFile("norm_2.fastq").toPath();
    FastqNormalizer.normalizePairedEnd(
        inputFile1.toString(),
        inputFile2.toString(),
        normOut1.toString(),
        normOut2.toString(),
        prefix,
        false,
        tempFolder.getRoot());

    List<FastqRecord> normRecords1 = readFastq(normOut1);
    List<FastqRecord> normRecords2 = readFastq(normOut2);

    assertTrue("Expected reads in pipeline output", pipelineOut1.size() > 0);
    assertEquals("Pair count mismatch /1", pipelineOut1.size(), normRecords1.size());
    assertEquals("Pair count mismatch /2", pipelineOut2.size(), normRecords2.size());

    // Build maps keyed by bases for comparison (ordering may differ slightly)
    Map<String, String> pipelineBasesToQual1 = new HashMap<>();
    for (FastqRecord r : pipelineOut1) {
      pipelineBasesToQual1.put(r.getReadString(), r.getBaseQualityString());
    }
    Map<String, String> normBasesToQual1 = new HashMap<>();
    for (FastqRecord r : normRecords1) {
      normBasesToQual1.put(r.getReadString(), r.getBaseQualityString());
    }

    assertEquals(
        "All base sequences should be present in both outputs",
        pipelineBasesToQual1.keySet(),
        normBasesToQual1.keySet());

    for (String bases : pipelineBasesToQual1.keySet()) {
      assertEquals(
          "Quality mismatch for bases: " + bases.substring(0, Math.min(20, bases.length())),
          pipelineBasesToQual1.get(bases),
          normBasesToQual1.get(bases));
    }
  }

  // ---- Single-end test ----

  /** Single-end: both paths should produce identical output. */
  @Test
  public void testSingleEnd() throws Exception {
    Path inputFile = tempFolder.newFile("single.fastq").toPath();
    Files.write(
        inputFile,
        ("@READ1\nACGTACGT\n+\nIIIIIIII\n"
                + "@READ2\nGGCCGGCC\n+\nJJJJJJJJ\n"
                + "@READ3\nTTAATTAA\n+\nKKKKKKKK\n")
            .getBytes(StandardCharsets.UTF_8));

    String prefix = "ERR888";

    // Run pipeline (single file → Fastq2Sam → Sam2Fastq)
    Path pipelineDir = tempFolder.newFolder("pipeline_se").toPath();

    Fastq2Sam.Params f2sParams = new Fastq2Sam.Params();
    f2sParams.tmp_root = pipelineDir.toString();
    f2sParams.sample_name = "SM-001";
    f2sParams.data_file = pipelineDir.resolve("output.bam").toString();
    f2sParams.compression = FileCompression.NONE.name();
    f2sParams.files = Collections.singletonList(inputFile.toString());

    Fastq2Sam f2s = new Fastq2Sam();
    f2s.create(f2sParams);

    Sam2Fastq.Params s2fParams = new Sam2Fastq.Params();
    s2fParams.samFile = new File(f2sParams.data_file);
    s2fParams.fastqBaseName = pipelineDir.resolve("result").toString();
    s2fParams.prefix = prefix;
    s2fParams.nofStreams = 3;

    Sam2Fastq s2f = new Sam2Fastq();
    s2f.create(s2fParams);

    // Pipeline single-end output goes to stream 0 (unpaired file)
    File pipelineOutFile = new File(s2fParams.fastqBaseName + ".fastq");
    List<FastqRecord> pipelineRecords = readFastq(pipelineOutFile.toPath());

    // Run normalizer
    Path normOut = tempFolder.newFile("norm_se.fastq").toPath();
    long normCount =
        FastqNormalizer.normalizeSingleEnd(inputFile.toString(), normOut.toString(), prefix, false);

    List<FastqRecord> normRecords = readFastq(normOut);

    assertEquals("Read count mismatch", pipelineRecords.size(), normRecords.size());
    assertEquals(3, normCount);

    // Sort by bases for stable comparison
    pipelineRecords.sort(Comparator.comparing(FastqRecord::getReadString));
    normRecords.sort(Comparator.comparing(FastqRecord::getReadString));

    for (int i = 0; i < pipelineRecords.size(); i++) {
      assertEquals(pipelineRecords.get(i).getReadString(), normRecords.get(i).getReadString());
      assertEquals(
          pipelineRecords.get(i).getBaseQualityString(), normRecords.get(i).getBaseQualityString());
    }
  }

  // ---- Helper methods ----

  /** Runs the Fastq2Sam → Sam2Fastq pipeline and collects output records. */
  private void runPipeline(
      Path inputFile1,
      Path inputFile2,
      String prefix,
      boolean convertUracil,
      List<FastqRecord> outRecords1,
      List<FastqRecord> outRecords2)
      throws Exception {

    Path pipelineDir = tempFolder.newFolder("pipeline").toPath();

    Fastq2Sam.Params f2sParams = new Fastq2Sam.Params();
    f2sParams.tmp_root = pipelineDir.toString();
    f2sParams.sample_name = "SM-001";
    f2sParams.data_file = pipelineDir.resolve("output.bam").toString();
    f2sParams.compression = FileCompression.NONE.name();
    f2sParams.convertUracil = convertUracil;
    f2sParams.files = Arrays.asList(inputFile1.toString(), inputFile2.toString());

    Fastq2Sam f2s = new Fastq2Sam();
    f2s.create(f2sParams);

    Sam2Fastq.Params s2fParams = new Sam2Fastq.Params();
    s2fParams.samFile = new File(f2sParams.data_file);
    s2fParams.fastqBaseName = pipelineDir.resolve("result").toString();
    s2fParams.prefix = prefix;
    s2fParams.nofStreams = 3;

    Sam2Fastq s2f = new Sam2Fastq();
    s2f.create(s2fParams);

    // Paired reads go to _1.fastq and _2.fastq
    File pipelineOutFile1 = new File(s2fParams.fastqBaseName + "_1.fastq");
    File pipelineOutFile2 = new File(s2fParams.fastqBaseName + "_2.fastq");

    // Unpaired reads go to .fastq (stream 0)
    File pipelineOutUnpaired = new File(s2fParams.fastqBaseName + ".fastq");

    if (pipelineOutFile1.exists()) {
      outRecords1.addAll(readFastq(pipelineOutFile1.toPath()));
    }
    if (pipelineOutFile2.exists()) {
      outRecords2.addAll(readFastq(pipelineOutFile2.toPath()));
    }
    // Append unpaired reads to outRecords1 (same as FastqNormalizer orphan behavior)
    if (pipelineOutUnpaired.exists()) {
      outRecords1.addAll(readFastq(pipelineOutUnpaired.toPath()));
    }
  }

  /**
   * Runs both paths for non-Casava paired reads and returns keyed results.
   *
   * @param reads Array of {name1, bases1, qual1, name2, bases2, qual2} per pair
   */
  private RoundTripResult runBothPaths(String[][] reads, String prefix, boolean convertUracil)
      throws Exception {

    Path inputFile1 = tempFolder.newFile("input_1.fastq").toPath();
    Path inputFile2 = tempFolder.newFile("input_2.fastq").toPath();

    StringBuilder sb1 = new StringBuilder();
    StringBuilder sb2 = new StringBuilder();
    for (String[] pair : reads) {
      sb1.append(pair[0]).append("\n").append(pair[1]).append("\n+\n").append(pair[2]).append("\n");
      sb2.append(pair[3]).append("\n").append(pair[4]).append("\n+\n").append(pair[5]).append("\n");
    }
    Files.write(inputFile1, sb1.toString().getBytes(StandardCharsets.UTF_8));
    Files.write(inputFile2, sb2.toString().getBytes(StandardCharsets.UTF_8));

    // Run pipeline
    List<FastqRecord> pipelineOut1 = new ArrayList<>();
    List<FastqRecord> pipelineOut2 = new ArrayList<>();
    runPipeline(inputFile1, inputFile2, prefix, convertUracil, pipelineOut1, pipelineOut2);

    // Run normalizer
    Path normOut1 = tempFolder.newFile("norm_1.fastq").toPath();
    Path normOut2 = tempFolder.newFile("norm_2.fastq").toPath();
    FastqNormalizer.normalizePairedEnd(
        inputFile1.toString(),
        inputFile2.toString(),
        normOut1.toString(),
        normOut2.toString(),
        prefix,
        convertUracil,
        tempFolder.getRoot());

    List<FastqRecord> normRecords1 = readFastq(normOut1);
    List<FastqRecord> normRecords2 = readFastq(normOut2);

    // Build keyed maps by extracting the read base name from the output.
    // Both paths output: PREFIX.N BASENAME/1 — extract BASENAME as key.
    Map<String, ReadPair> pipelinePairs = buildPairMap(pipelineOut1, pipelineOut2);
    Map<String, ReadPair> normPairs = buildPairMap(normRecords1, normRecords2);

    return new RoundTripResult(pipelinePairs, normPairs);
  }

  /**
   * Builds a map of read key → ReadPair from paired output files. Assumes read names are in format
   * "PREFIX.N BASENAME/1" or "PREFIX.N BASENAME/2".
   */
  private Map<String, ReadPair> buildPairMap(
      List<FastqRecord> records1, List<FastqRecord> records2) {
    assertEquals(
        "File 1 and file 2 should have same number of records", records1.size(), records2.size());

    // Sort both by the base name portion for stable matching
    Comparator<FastqRecord> byBaseName =
        Comparator.comparing(r -> extractBaseName(r.getReadName()));
    records1.sort(byBaseName);
    records2.sort(byBaseName);

    Map<String, ReadPair> pairs = new LinkedHashMap<>();
    for (int i = 0; i < records1.size(); i++) {
      String key = extractBaseName(records1.get(i).getReadName());
      ReadPair pair = new ReadPair();
      pair.bases1 = records1.get(i).getReadString();
      pair.qual1 = records1.get(i).getBaseQualityString();
      pair.bases2 = records2.get(i).getReadString();
      pair.qual2 = records2.get(i).getBaseQualityString();
      pairs.put(key, pair);
    }
    return pairs;
  }

  /**
   * Extracts the base name from an output read name. Input formats:
   *
   * <ul>
   *   <li>"PREFIX.N BASENAME/1" → "BASENAME"
   *   <li>"PREFIX.N BASENAME" → "BASENAME" (orphan)
   * </ul>
   */
  private String extractBaseName(String readName) {
    // Strip "PREFIX.N " part
    int spaceIdx = readName.indexOf(' ');
    String afterPrefix = spaceIdx >= 0 ? readName.substring(spaceIdx + 1) : readName;
    // Strip /1 or /2 suffix
    if (afterPrefix.endsWith("/1") || afterPrefix.endsWith("/2")) {
      afterPrefix = afterPrefix.substring(0, afterPrefix.length() - 2);
    }
    return afterPrefix;
  }

  private List<FastqRecord> readFastq(Path path) throws IOException {
    List<FastqRecord> records = new ArrayList<>();
    try (FastqReader reader =
        new FastqReader(
            null,
            new BufferedReader(
                new InputStreamReader(Utils.openFastqInputStream(path), StandardCharsets.UTF_8)))) {
      for (FastqRecord record : reader) {
        records.add(record);
      }
    }
    return records;
  }

  private static class ReadPair {
    String bases1, qual1;
    String bases2, qual2;
  }

  private static class RoundTripResult {
    final Map<String, ReadPair> pipelinePairs;
    final Map<String, ReadPair> normalizerPairs;

    RoundTripResult(Map<String, ReadPair> pipelinePairs, Map<String, ReadPair> normalizerPairs) {
      this.pipelinePairs = pipelinePairs;
      this.normalizerPairs = normalizerPairs;
    }
  }
}
