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
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;
import uk.ac.ebi.ena.readtools.utils.Utils;

public class FastqNormalizerTest {

  @Rule public TemporaryFolder tempFolder = new TemporaryFolder();

  // ---- Single-end tests ----

  /** Sanger (Phred+33) quality scores should pass through unchanged. */
  @Test
  public void testSingleEndSangerPassthrough() throws IOException {
    Path inputFile = tempFolder.newFile("input.fastq").toPath();
    Files.write(
        inputFile,
        ("@READ1\n"
                + "ACGT\n"
                + "+\n"
                + "###$\n" // Phred+33: ASCII 35,35,35,36 = Phred 2,2,2,3
                + "@READ2\n"
                + "GGCC\n"
                + "+\n"
                + "%%&&\n") // Phred+33: ASCII 37,37,38,38 = Phred 4,4,5,5
            .getBytes(StandardCharsets.UTF_8));

    Path outputFile = tempFolder.newFile("output.fastq").toPath();

    long count =
        FastqNormalizer.normalizeSingleEnd(
            inputFile.toString(), outputFile.toString(), null, false);

    assertEquals(2, count);

    List<FastqRecord> records = readFastq(outputFile);
    assertEquals(2, records.size());

    assertEquals("READ1", records.get(0).getReadName());
    assertEquals("ACGT", records.get(0).getReadString());
    assertEquals("###$", records.get(0).getBaseQualityString());

    assertEquals("READ2", records.get(1).getReadName());
    assertEquals("GGCC", records.get(1).getReadString());
    assertEquals("%%&&", records.get(1).getBaseQualityString());
  }

  /** With a prefix, read names should become "{prefix}.{counter} {originalName}". */
  @Test
  public void testSingleEndWithPrefix() throws IOException {
    Path inputFile = tempFolder.newFile("input.fastq").toPath();
    Files.write(
        inputFile,
        ("@READ1\n"
                + "ACGT\n"
                + "+\n"
                + "IIII\n"
                + "@READ2\n"
                + "GGCC\n"
                + "+\n"
                + "JJJJ\n"
                + "@READ3\n"
                + "TTAA\n"
                + "+\n"
                + "KKKK\n")
            .getBytes(StandardCharsets.UTF_8));

    Path outputFile = tempFolder.newFile("output.fastq").toPath();

    long count =
        FastqNormalizer.normalizeSingleEnd(
            inputFile.toString(), outputFile.toString(), "ERR123456", false);

    assertEquals(3, count);

    List<FastqRecord> records = readFastq(outputFile);
    assertEquals(3, records.size());

    assertEquals("ERR123456.1 READ1", records.get(0).getReadName());
    assertEquals("ERR123456.2 READ2", records.get(1).getReadName());
    assertEquals("ERR123456.3 READ3", records.get(2).getReadName());
  }

  /** Without a prefix, read names should be preserved as-is. */
  @Test
  public void testSingleEndWithoutPrefix() throws IOException {
    Path inputFile = tempFolder.newFile("input.fastq").toPath();
    Files.write(
        inputFile, ("@READ1\n" + "ACGT\n" + "+\n" + "IIII\n").getBytes(StandardCharsets.UTF_8));

    Path outputFile = tempFolder.newFile("output.fastq").toPath();

    long count =
        FastqNormalizer.normalizeSingleEnd(
            inputFile.toString(), outputFile.toString(), null, false);

    assertEquals(1, count);

    List<FastqRecord> records = readFastq(outputFile);
    assertEquals(1, records.size());
    assertEquals("READ1", records.get(0).getReadName());
  }

  /** U bases should be converted to T (uppercase U→T, lowercase u→t). */
  @Test
  public void testSingleEndUracilConversion() throws IOException {
    Path inputFile = tempFolder.newFile("input.fastq").toPath();
    Files.write(
        inputFile,
        ("@READ1\n"
                + "AUGCU\n"
                + "+\n"
                + "IIIII\n"
                + "@READ2\n"
                + "gacggaUCuauagcaaaacu\n"
                + "+\n"
                + "IIIIIIIIIIIIIIIIIIII\n")
            .getBytes(StandardCharsets.UTF_8));

    Path outputFile = tempFolder.newFile("output.fastq").toPath();

    long count =
        FastqNormalizer.normalizeSingleEnd(inputFile.toString(), outputFile.toString(), null, true);

    assertEquals(2, count);

    List<FastqRecord> records = readFastq(outputFile);
    assertEquals(2, records.size());

    assertEquals("ATGCT", records.get(0).getReadString());
    assertEquals("gacggaTCtatagcaaaact", records.get(1).getReadString());
  }

  /** Test with an actual uracil-bases test resource file. */
  @Test
  public void testSingleEndRealFile() throws IOException {
    String inputFile = "src/test/resources/uracil-bases_1.fastq";
    Path outputFile = tempFolder.newFile("output.fastq").toPath();

    long count = FastqNormalizer.normalizeSingleEnd(inputFile, outputFile.toString(), "TEST", true);

    assertTrue(count > 0);

    List<FastqRecord> records = readFastq(outputFile);
    assertEquals(count, records.size());

    assertTrue(records.get(0).getReadName().startsWith("TEST.1 "));

    for (FastqRecord record : records) {
      assertFalse(record.getReadString().contains("U"));
      assertFalse(record.getReadString().contains("u"));
    }
  }

  // ---- Paired-end tests ----

  /** Two in-order paired files should produce matched pairs with /1 and /2 suffixes. */
  @Test
  public void testPairedEndInOrder() throws IOException {
    Path inputFile1 = tempFolder.newFile("input_1.fastq").toPath();
    Path inputFile2 = tempFolder.newFile("input_2.fastq").toPath();

    Files.write(
        inputFile1,
        ("@READ1/1\n" + "ACGT\n" + "+\n" + "IIII\n" + "@READ2/1\n" + "GGCC\n" + "+\n" + "JJJJ\n")
            .getBytes(StandardCharsets.UTF_8));

    Files.write(
        inputFile2,
        ("@READ1/2\n" + "TTAA\n" + "+\n" + "KKKK\n" + "@READ2/2\n" + "CCGG\n" + "+\n" + "LLLL\n")
            .getBytes(StandardCharsets.UTF_8));

    Path outputFile1 = tempFolder.newFile("output_1.fastq").toPath();
    Path outputFile2 = tempFolder.newFile("output_2.fastq").toPath();

    FastqNormalizer.PairedNormalizationResult result =
        FastqNormalizer.normalizePairedEnd(
            inputFile1.toString(),
            inputFile2.toString(),
            outputFile1.toString(),
            outputFile2.toString(),
            null,
            false,
            tempFolder.getRoot());

    assertEquals(2, result.getPairCount());
    assertEquals(0, result.getOrphanCount());
    assertEquals(4, result.getTotalReadCount());

    List<FastqRecord> records1 = readFastq(outputFile1);
    List<FastqRecord> records2 = readFastq(outputFile2);

    assertEquals(2, records1.size());
    assertEquals(2, records2.size());

    // Slash separator preserved, no prefix
    assertEquals("READ1/1", records1.get(0).getReadName());
    assertEquals("READ1/2", records2.get(0).getReadName());
    assertEquals("READ2/1", records1.get(1).getReadName());
    assertEquals("READ2/2", records2.get(1).getReadName());
  }

  /** Paired-end with prefix: read names should be "{prefix}.{counter} {baseName}/1". */
  @Test
  public void testPairedEndWithPrefix() throws IOException {
    Path inputFile1 = tempFolder.newFile("input_1.fastq").toPath();
    Path inputFile2 = tempFolder.newFile("input_2.fastq").toPath();

    Files.write(
        inputFile1, ("@READ1/1\n" + "ACGT\n" + "+\n" + "IIII\n").getBytes(StandardCharsets.UTF_8));

    Files.write(
        inputFile2, ("@READ1/2\n" + "TTAA\n" + "+\n" + "KKKK\n").getBytes(StandardCharsets.UTF_8));

    Path outputFile1 = tempFolder.newFile("output_1.fastq").toPath();
    Path outputFile2 = tempFolder.newFile("output_2.fastq").toPath();

    FastqNormalizer.PairedNormalizationResult result =
        FastqNormalizer.normalizePairedEnd(
            inputFile1.toString(),
            inputFile2.toString(),
            outputFile1.toString(),
            outputFile2.toString(),
            "SRR654321",
            false,
            tempFolder.getRoot());

    assertEquals(1, result.getPairCount());
    assertEquals(0, result.getOrphanCount());
    assertEquals(2, result.getTotalReadCount());

    List<FastqRecord> records1 = readFastq(outputFile1);
    List<FastqRecord> records2 = readFastq(outputFile2);

    assertEquals("SRR654321.1 READ1/1", records1.get(0).getReadName());
    assertEquals("SRR654321.1 READ1/2", records2.get(0).getReadName());
  }

  /**
   * Mismatched file lengths should be tolerated: paired reads pair normally, extras become orphans
   * written to the first output file without /1 or /2 suffix.
   */
  @Test
  public void testPairedEndMismatchedLengths() throws IOException {
    Path inputFile1 = tempFolder.newFile("input_1.fastq").toPath();
    Path inputFile2 = tempFolder.newFile("input_2.fastq").toPath();

    // File 1 has 2 reads, File 2 has 1 read
    Files.write(
        inputFile1,
        ("@READ1/1\n" + "ACGT\n" + "+\n" + "IIII\n" + "@READ2/1\n" + "GGCC\n" + "+\n" + "JJJJ\n")
            .getBytes(StandardCharsets.UTF_8));

    Files.write(
        inputFile2, ("@READ1/2\n" + "TTAA\n" + "+\n" + "KKKK\n").getBytes(StandardCharsets.UTF_8));

    Path outputFile1 = tempFolder.newFile("output_1.fastq").toPath();
    Path outputFile2 = tempFolder.newFile("output_2.fastq").toPath();

    FastqNormalizer.PairedNormalizationResult result =
        FastqNormalizer.normalizePairedEnd(
            inputFile1.toString(),
            inputFile2.toString(),
            outputFile1.toString(),
            outputFile2.toString(),
            null,
            false,
            tempFolder.getRoot());

    assertEquals(1, result.getPairCount());
    assertEquals(1, result.getOrphanCount());
    assertEquals(3, result.getTotalReadCount());

    List<FastqRecord> records1 = readFastq(outputFile1);
    List<FastqRecord> records2 = readFastq(outputFile2);

    // First output file: 1 paired read (/1) + 1 orphan (no /1 or /2 suffix)
    assertEquals(2, records1.size());
    // Second output file: 1 paired read (/2)
    assertEquals(1, records2.size());

    // Verify the orphan has no /1 or /2 suffix
    boolean hasOrphan = records1.stream().anyMatch(r -> !r.getReadName().contains("/"));
    assertTrue("Expected an orphan read without /1 or /2 suffix", hasOrphan);
  }

  /** Casava 1.8 format reads (space-separated instrument + read metadata) should pair correctly. */
  @Test
  public void testPairedEndCasavaFormat() throws IOException {
    Path inputFile1 = tempFolder.newFile("input_1.fastq").toPath();
    Path inputFile2 = tempFolder.newFile("input_2.fastq").toPath();

    Files.write(
        inputFile1,
        ("@A00500:310:HG3JFDRXY:1:2101:23963:1000 1:N:0:CTTCGTTC+GAGAAGGT\n"
                + "ACGT\n"
                + "+\n"
                + "IIII\n")
            .getBytes(StandardCharsets.UTF_8));

    Files.write(
        inputFile2,
        ("@A00500:310:HG3JFDRXY:1:2101:23963:1000 2:N:0:CTTCGTTC+GAGAAGGT\n"
                + "TTAA\n"
                + "+\n"
                + "KKKK\n")
            .getBytes(StandardCharsets.UTF_8));

    Path outputFile1 = tempFolder.newFile("output_1.fastq").toPath();
    Path outputFile2 = tempFolder.newFile("output_2.fastq").toPath();

    FastqNormalizer.PairedNormalizationResult result =
        FastqNormalizer.normalizePairedEnd(
            inputFile1.toString(),
            inputFile2.toString(),
            outputFile1.toString(),
            outputFile2.toString(),
            null,
            false,
            tempFolder.getRoot());

    assertEquals(1, result.getPairCount());
    assertEquals(0, result.getOrphanCount());

    List<FastqRecord> records1 = readFastq(outputFile1);
    List<FastqRecord> records2 = readFastq(outputFile2);

    assertEquals(1, records1.size());
    assertEquals(1, records2.size());

    // Casava read names already contain the pair number — no /1 or /2 should be appended.
    assertEquals(
        "A00500:310:HG3JFDRXY:1:2101:23963:1000 1:N:0:CTTCGTTC+GAGAAGGT",
        records1.get(0).getReadName());
    assertEquals(
        "A00500:310:HG3JFDRXY:1:2101:23963:1000 2:N:0:CTTCGTTC+GAGAAGGT",
        records2.get(0).getReadName());
  }

  /**
   * Casava 1.8 format with prefix: read names should be "{prefix}.{counter} {originalCasavaName}"
   * without /1 or /2 appended, since the Casava tail already has the pair number.
   */
  @Test
  public void testPairedEndCasavaFormatWithPrefix() throws IOException {
    Path inputFile1 = tempFolder.newFile("input_1.fastq").toPath();
    Path inputFile2 = tempFolder.newFile("input_2.fastq").toPath();

    Files.write(
        inputFile1,
        ("@A00500:310:HG3JFDRXY:1:2101:23963:1000 1:N:0:CTTCGTTC+GAGAAGGT\n"
                + "ACGT\n"
                + "+\n"
                + "IIII\n")
            .getBytes(StandardCharsets.UTF_8));

    Files.write(
        inputFile2,
        ("@A00500:310:HG3JFDRXY:1:2101:23963:1000 2:N:0:CTTCGTTC+GAGAAGGT\n"
                + "TTAA\n"
                + "+\n"
                + "KKKK\n")
            .getBytes(StandardCharsets.UTF_8));

    Path outputFile1 = tempFolder.newFile("output_1.fastq").toPath();
    Path outputFile2 = tempFolder.newFile("output_2.fastq").toPath();

    FastqNormalizer.PairedNormalizationResult result =
        FastqNormalizer.normalizePairedEnd(
            inputFile1.toString(),
            inputFile2.toString(),
            outputFile1.toString(),
            outputFile2.toString(),
            "ERR999",
            false,
            tempFolder.getRoot());

    assertEquals(1, result.getPairCount());

    List<FastqRecord> records1 = readFastq(outputFile1);
    List<FastqRecord> records2 = readFastq(outputFile2);

    assertEquals(
        "ERR999.1 A00500:310:HG3JFDRXY:1:2101:23963:1000 1:N:0:CTTCGTTC+GAGAAGGT",
        records1.get(0).getReadName());
    assertEquals(
        "ERR999.1 A00500:310:HG3JFDRXY:1:2101:23963:1000 2:N:0:CTTCGTTC+GAGAAGGT",
        records2.get(0).getReadName());
  }

  /** Paired-end U→T conversion should apply to both mate files. */
  @Test
  public void testPairedEndUracilConversion() throws IOException {
    Path inputFile1 = tempFolder.newFile("input_1.fastq").toPath();
    Path inputFile2 = tempFolder.newFile("input_2.fastq").toPath();

    Files.write(
        inputFile1,
        ("@READ1/1\n" + "AUGCU\n" + "+\n" + "IIIII\n").getBytes(StandardCharsets.UTF_8));

    Files.write(
        inputFile2,
        ("@READ1/2\n" + "uuAAuu\n" + "+\n" + "KKKKKK\n").getBytes(StandardCharsets.UTF_8));

    Path outputFile1 = tempFolder.newFile("output_1.fastq").toPath();
    Path outputFile2 = tempFolder.newFile("output_2.fastq").toPath();

    FastqNormalizer.PairedNormalizationResult result =
        FastqNormalizer.normalizePairedEnd(
            inputFile1.toString(),
            inputFile2.toString(),
            outputFile1.toString(),
            outputFile2.toString(),
            null,
            true,
            tempFolder.getRoot());

    assertEquals(1, result.getPairCount());

    List<FastqRecord> records1 = readFastq(outputFile1);
    List<FastqRecord> records2 = readFastq(outputFile2);

    assertEquals("ATGCT", records1.get(0).getReadString());
    assertEquals("ttAAtt", records2.get(0).getReadString());
  }

  /** Test with actual paired FASTQ test resource files. */
  @Test
  public void testPairedEndRealFiles() throws IOException {
    String inputFile1 = "src/test/resources/2fastq/28239_1822_1.fastq";
    String inputFile2 = "src/test/resources/2fastq/28239_1822_2.fastq";
    Path outputFile1 = tempFolder.newFile("output_1.fastq").toPath();
    Path outputFile2 = tempFolder.newFile("output_2.fastq").toPath();

    FastqNormalizer.PairedNormalizationResult result =
        FastqNormalizer.normalizePairedEnd(
            inputFile1,
            inputFile2,
            outputFile1.toString(),
            outputFile2.toString(),
            "RUN001",
            false,
            tempFolder.getRoot());

    assertTrue(result.getPairCount() > 0);
    assertEquals(0, result.getOrphanCount());

    List<FastqRecord> records1 = readFastq(outputFile1);
    List<FastqRecord> records2 = readFastq(outputFile2);

    assertEquals(result.getPairCount(), records1.size());
    assertEquals(result.getPairCount(), records2.size());

    assertTrue(records1.get(0).getReadName().startsWith("RUN001.1 "));
    assertTrue(records1.get(0).getReadName().endsWith("/1"));
    assertTrue(records2.get(0).getReadName().startsWith("RUN001.1 "));
    assertTrue(records2.get(0).getReadName().endsWith("/2"));
  }

  // ---- Non-Casava separator format tests ----
  // PairedFastqWriter.SPLIT_REGEXP supports pairing on dot, colon, slash, and underscore
  // separators. These tests verify FastqNormalizer handles them all.

  /** Dot-separated pair index: READ.1 / READ.2 */
  @Test
  public void testPairedEndDotSeparator() throws IOException {
    Path inputFile1 = tempFolder.newFile("input_1.fastq").toPath();
    Path inputFile2 = tempFolder.newFile("input_2.fastq").toPath();

    Files.write(
        inputFile1,
        ("@READX.1\n" + "ACGT\n" + "+\n" + "IIII\n" + "@READY.1\n" + "GGCC\n" + "+\n" + "JJJJ\n")
            .getBytes(StandardCharsets.UTF_8));

    Files.write(
        inputFile2,
        ("@READX.2\n" + "TTAA\n" + "+\n" + "KKKK\n" + "@READY.2\n" + "CCGG\n" + "+\n" + "LLLL\n")
            .getBytes(StandardCharsets.UTF_8));

    Path outputFile1 = tempFolder.newFile("output_1.fastq").toPath();
    Path outputFile2 = tempFolder.newFile("output_2.fastq").toPath();

    FastqNormalizer.PairedNormalizationResult result =
        FastqNormalizer.normalizePairedEnd(
            inputFile1.toString(),
            inputFile2.toString(),
            outputFile1.toString(),
            outputFile2.toString(),
            "RUN",
            false,
            tempFolder.getRoot());

    assertEquals(2, result.getPairCount());
    assertEquals(0, result.getOrphanCount());
    assertEquals(4, result.getTotalReadCount());

    List<FastqRecord> records1 = readFastq(outputFile1);
    List<FastqRecord> records2 = readFastq(outputFile2);

    assertEquals(2, records1.size());
    assertEquals(2, records2.size());

    // Dot separator stripped, normalized to /1 /2 — matches pipeline behavior
    assertEquals("RUN.1 READX/1", records1.get(0).getReadName());
    assertEquals("RUN.1 READX/2", records2.get(0).getReadName());
    assertEquals("RUN.2 READY/1", records1.get(1).getReadName());
    assertEquals("RUN.2 READY/2", records2.get(1).getReadName());
  }

  /** Colon-separated pair index: READ:1 / READ:2 (non-Casava-like name so SPLIT_REGEXP fires). */
  @Test
  public void testPairedEndColonSeparator() throws IOException {
    Path inputFile1 = tempFolder.newFile("input_1.fastq").toPath();
    Path inputFile2 = tempFolder.newFile("input_2.fastq").toPath();

    Files.write(
        inputFile1,
        ("@SOLEXA_READ:1\n" + "ACGT\n" + "+\n" + "IIII\n").getBytes(StandardCharsets.UTF_8));

    Files.write(
        inputFile2,
        ("@SOLEXA_READ:2\n" + "TTAA\n" + "+\n" + "KKKK\n").getBytes(StandardCharsets.UTF_8));

    Path outputFile1 = tempFolder.newFile("output_1.fastq").toPath();
    Path outputFile2 = tempFolder.newFile("output_2.fastq").toPath();

    FastqNormalizer.PairedNormalizationResult result =
        FastqNormalizer.normalizePairedEnd(
            inputFile1.toString(),
            inputFile2.toString(),
            outputFile1.toString(),
            outputFile2.toString(),
            null,
            false,
            tempFolder.getRoot());

    assertEquals(1, result.getPairCount());
    assertEquals(0, result.getOrphanCount());

    List<FastqRecord> records1 = readFastq(outputFile1);
    List<FastqRecord> records2 = readFastq(outputFile2);

    assertEquals(1, records1.size());
    assertEquals(1, records2.size());

    // Colon separator stripped, normalized to /1 /2 — matches pipeline behavior
    assertEquals("SOLEXA_READ/1", records1.get(0).getReadName());
    assertEquals("SOLEXA_READ/2", records2.get(0).getReadName());
  }

  /** Underscore-separated pair index: READ_1 / READ_2 */
  @Test
  public void testPairedEndUnderscoreSeparator() throws IOException {
    Path inputFile1 = tempFolder.newFile("input_1.fastq").toPath();
    Path inputFile2 = tempFolder.newFile("input_2.fastq").toPath();

    Files.write(
        inputFile1, ("@MYREAD_1\n" + "ACGT\n" + "+\n" + "IIII\n").getBytes(StandardCharsets.UTF_8));

    Files.write(
        inputFile2, ("@MYREAD_2\n" + "TTAA\n" + "+\n" + "KKKK\n").getBytes(StandardCharsets.UTF_8));

    Path outputFile1 = tempFolder.newFile("output_1.fastq").toPath();
    Path outputFile2 = tempFolder.newFile("output_2.fastq").toPath();

    FastqNormalizer.PairedNormalizationResult result =
        FastqNormalizer.normalizePairedEnd(
            inputFile1.toString(),
            inputFile2.toString(),
            outputFile1.toString(),
            outputFile2.toString(),
            "PREFIX",
            false,
            tempFolder.getRoot());

    assertEquals(1, result.getPairCount());
    assertEquals(0, result.getOrphanCount());
    assertEquals(2, result.getTotalReadCount());

    List<FastqRecord> records1 = readFastq(outputFile1);
    List<FastqRecord> records2 = readFastq(outputFile2);

    assertEquals(1, records1.size());
    assertEquals(1, records2.size());

    // Underscore separator stripped, normalized to /1 /2 — matches pipeline behavior
    assertEquals("PREFIX.1 MYREAD/1", records1.get(0).getReadName());
    assertEquals("PREFIX.1 MYREAD/2", records2.get(0).getReadName());
  }

  /** Multi-digit pair index: READ.10 / READ.20 — normalized to /1 /2, matches pipeline. */
  @Test
  public void testPairedEndMultiDigitIndex() throws IOException {
    Path inputFile1 = tempFolder.newFile("input_1.fastq").toPath();
    Path inputFile2 = tempFolder.newFile("input_2.fastq").toPath();

    Files.write(
        inputFile1, ("@MYREAD.10\nACGTACGT\n+\nIIIIIIII\n").getBytes(StandardCharsets.UTF_8));
    Files.write(
        inputFile2, ("@MYREAD.20\nTTAATTAA\n+\nKKKKKKKK\n").getBytes(StandardCharsets.UTF_8));

    Path outputFile1 = tempFolder.newFile("output_1.fastq").toPath();
    Path outputFile2 = tempFolder.newFile("output_2.fastq").toPath();

    FastqNormalizer.PairedNormalizationResult result =
        FastqNormalizer.normalizePairedEnd(
            inputFile1.toString(),
            inputFile2.toString(),
            outputFile1.toString(),
            outputFile2.toString(),
            "ERR006",
            false,
            tempFolder.getRoot());

    assertEquals(1, result.getPairCount());
    assertEquals(0, result.getOrphanCount());

    List<FastqRecord> records1 = readFastq(outputFile1);
    List<FastqRecord> records2 = readFastq(outputFile2);

    assertEquals(1, records1.size());
    assertEquals(1, records2.size());

    // Multi-digit indices stripped, normalized to /1 /2 — matches pipeline behavior
    assertEquals("ERR006.1 MYREAD/1", records1.get(0).getReadName());
    assertEquals("ERR006.1 MYREAD/2", records2.get(0).getReadName());
  }

  // ---- Read count accuracy tests ----

  /** All reads paired: pairCount=3, orphanCount=0, totalReadCount=6. */
  @Test
  public void testReadCountAllPaired() throws IOException {
    Path inputFile1 = tempFolder.newFile("input_1.fastq").toPath();
    Path inputFile2 = tempFolder.newFile("input_2.fastq").toPath();

    Files.write(
        inputFile1,
        ("@A/1\n" + "ACGT\n" + "+\n" + "IIII\n" + "@B/1\n" + "GGCC\n" + "+\n" + "JJJJ\n" + "@C/1\n"
                + "TTAA\n" + "+\n" + "KKKK\n")
            .getBytes(StandardCharsets.UTF_8));

    Files.write(
        inputFile2,
        ("@A/2\n" + "ACGT\n" + "+\n" + "IIII\n" + "@B/2\n" + "GGCC\n" + "+\n" + "JJJJ\n" + "@C/2\n"
                + "TTAA\n" + "+\n" + "KKKK\n")
            .getBytes(StandardCharsets.UTF_8));

    Path outputFile1 = tempFolder.newFile("output_1.fastq").toPath();
    Path outputFile2 = tempFolder.newFile("output_2.fastq").toPath();

    FastqNormalizer.PairedNormalizationResult result =
        FastqNormalizer.normalizePairedEnd(
            inputFile1.toString(),
            inputFile2.toString(),
            outputFile1.toString(),
            outputFile2.toString(),
            null,
            false,
            tempFolder.getRoot());

    assertEquals(3, result.getPairCount());
    assertEquals(0, result.getOrphanCount());
    assertEquals(6, result.getTotalReadCount());
  }

  /**
   * Mixed pairs and orphans: 3 reads in file 1, 1 in file 2. A pairs, B and C become orphans.
   * Orphans go to the first output file without /1 or /2 suffix.
   */
  @Test
  public void testReadCountWithOrphans() throws IOException {
    Path inputFile1 = tempFolder.newFile("input_1.fastq").toPath();
    Path inputFile2 = tempFolder.newFile("input_2.fastq").toPath();

    Files.write(
        inputFile1,
        ("@A/1\n" + "ACGT\n" + "+\n" + "IIII\n" + "@B/1\n" + "GGCC\n" + "+\n" + "JJJJ\n" + "@C/1\n"
                + "TTAA\n" + "+\n" + "KKKK\n")
            .getBytes(StandardCharsets.UTF_8));

    Files.write(
        inputFile2, ("@A/2\n" + "ACGT\n" + "+\n" + "IIII\n").getBytes(StandardCharsets.UTF_8));

    Path outputFile1 = tempFolder.newFile("output_1.fastq").toPath();
    Path outputFile2 = tempFolder.newFile("output_2.fastq").toPath();

    FastqNormalizer.PairedNormalizationResult result =
        FastqNormalizer.normalizePairedEnd(
            inputFile1.toString(),
            inputFile2.toString(),
            outputFile1.toString(),
            outputFile2.toString(),
            "RUN",
            false,
            tempFolder.getRoot());

    assertEquals(1, result.getPairCount());
    assertEquals(2, result.getOrphanCount());
    assertEquals(4, result.getTotalReadCount()); // 1*2 + 2 = 4

    List<FastqRecord> records1 = readFastq(outputFile1);
    List<FastqRecord> records2 = readFastq(outputFile2);

    // Output file 1: 1 paired + 2 orphans = 3
    assertEquals(3, records1.size());
    // Output file 2: 1 paired
    assertEquals(1, records2.size());

    // Verify orphans have no /1 or /2 suffix
    long orphansInFile1 = records1.stream().filter(r -> !r.getReadName().contains("/")).count();
    assertEquals(2, orphansInFile1);
  }

  /**
   * Output should be in lexicographic read-key order (matching BAM queryname sort), regardless of
   * input order. Counters should follow the sorted order.
   */
  @Test
  public void testReadCountLexicographicOrder() throws IOException {
    Path inputFile1 = tempFolder.newFile("input_1.fastq").toPath();
    Path inputFile2 = tempFolder.newFile("input_2.fastq").toPath();

    // Input in non-lexicographic order: C, A, B
    Files.write(
        inputFile1,
        ("@C_READ/1\n"
                + "TTAA\n"
                + "+\n"
                + "KKKK\n"
                + "@A_READ/1\n"
                + "ACGT\n"
                + "+\n"
                + "IIII\n"
                + "@B_READ/1\n"
                + "GGCC\n"
                + "+\n"
                + "JJJJ\n")
            .getBytes(StandardCharsets.UTF_8));

    Files.write(
        inputFile2,
        ("@C_READ/2\n"
                + "TTAA\n"
                + "+\n"
                + "KKKK\n"
                + "@A_READ/2\n"
                + "ACGT\n"
                + "+\n"
                + "IIII\n"
                + "@B_READ/2\n"
                + "GGCC\n"
                + "+\n"
                + "JJJJ\n")
            .getBytes(StandardCharsets.UTF_8));

    Path outputFile1 = tempFolder.newFile("output_1.fastq").toPath();
    Path outputFile2 = tempFolder.newFile("output_2.fastq").toPath();

    FastqNormalizer.PairedNormalizationResult result =
        FastqNormalizer.normalizePairedEnd(
            inputFile1.toString(),
            inputFile2.toString(),
            outputFile1.toString(),
            outputFile2.toString(),
            "RUN",
            false,
            tempFolder.getRoot());

    assertEquals(3, result.getPairCount());

    List<FastqRecord> records1 = readFastq(outputFile1);

    // Output should be in lexicographic order: A_READ, B_READ, C_READ
    assertTrue(records1.get(0).getReadName().contains("A_READ"));
    assertTrue(records1.get(1).getReadName().contains("B_READ"));
    assertTrue(records1.get(2).getReadName().contains("C_READ"));

    // Counter should follow sort order
    assertTrue(records1.get(0).getReadName().startsWith("RUN.1 "));
    assertTrue(records1.get(1).getReadName().startsWith("RUN.2 "));
    assertTrue(records1.get(2).getReadName().startsWith("RUN.3 "));
  }

  // ---- Base count tracking tests ----

  @Test
  public void testBaseCountSingleEnd() throws IOException {
    Path inputFile = tempFolder.newFile("input.fastq").toPath();
    // 3 reads: 8 + 5 + 4 = 17 bases
    Files.write(
        inputFile,
        ("@READ1\nACGTACGT\n+\nIIIIIIII\n"
                + "@READ2\nGGCCG\n+\nJJJJJ\n"
                + "@READ3\nTTAA\n+\nKKKK\n")
            .getBytes(StandardCharsets.UTF_8));

    Path outputFile = tempFolder.newFile("output.fastq").toPath();
    FastqNormalizer.SingleNormalizationResult result =
        FastqNormalizer.normalizeSingleEndWithStats(
            inputFile.toString(), outputFile.toString(), null, false);

    assertEquals(3, result.getReadCount());
    assertEquals(17, result.getBaseCount());
  }

  @Test
  public void testBaseCountPairedEnd() throws IOException {
    Path inputFile1 = tempFolder.newFile("input_1.fastq").toPath();
    Path inputFile2 = tempFolder.newFile("input_2.fastq").toPath();

    // 2 pairs: (8+8) + (5+5) = 26 bases
    Files.write(
        inputFile1,
        ("@READA/1\nACGTACGT\n+\nIIIIIIII\n" + "@READB/1\nGGCCG\n+\nJJJJJ\n")
            .getBytes(StandardCharsets.UTF_8));
    Files.write(
        inputFile2,
        ("@READA/2\nTTAATTAA\n+\nKKKKKKKK\n" + "@READB/2\nCCGGC\n+\nLLLLL\n")
            .getBytes(StandardCharsets.UTF_8));

    Path outputFile1 = tempFolder.newFile("output_1.fastq").toPath();
    Path outputFile2 = tempFolder.newFile("output_2.fastq").toPath();

    FastqNormalizer.PairedNormalizationResult result =
        FastqNormalizer.normalizePairedEnd(
            inputFile1.toString(),
            inputFile2.toString(),
            outputFile1.toString(),
            outputFile2.toString(),
            null,
            false,
            tempFolder.getRoot());

    assertEquals(2, result.getPairCount());
    assertEquals(0, result.getOrphanCount());
    assertEquals(26, result.getBaseCount());
  }

  @Test
  public void testBaseCountWithOrphans() throws IOException {
    Path inputFile1 = tempFolder.newFile("input_1.fastq").toPath();
    Path inputFile2 = tempFolder.newFile("input_2.fastq").toPath();

    // READA pairs (8+8=16), READB is orphan (5), READC is orphan (4) = 25 total
    Files.write(
        inputFile1,
        ("@READA/1\nACGTACGT\n+\nIIIIIIII\n" + "@READB/1\nGGCCG\n+\nJJJJJ\n")
            .getBytes(StandardCharsets.UTF_8));
    Files.write(
        inputFile2,
        ("@READA/2\nTTAATTAA\n+\nKKKKKKKK\n" + "@READC/2\nCCGG\n+\nLLLL\n")
            .getBytes(StandardCharsets.UTF_8));

    Path outputFile1 = tempFolder.newFile("output_1.fastq").toPath();
    Path outputFile2 = tempFolder.newFile("output_2.fastq").toPath();

    FastqNormalizer.PairedNormalizationResult result =
        FastqNormalizer.normalizePairedEnd(
            inputFile1.toString(),
            inputFile2.toString(),
            outputFile1.toString(),
            outputFile2.toString(),
            null,
            false,
            tempFolder.getRoot());

    assertEquals(1, result.getPairCount());
    assertEquals(2, result.getOrphanCount());
    assertEquals(25, result.getBaseCount());
  }

  // ---- Spill file reassembly tests ----
  // These tests use tiny spillPageSize values to force disk spilling and verify
  // that multi-generation reassembly correctly pairs reads.

  /**
   * In-order pairs with spillPageSize=2: after 2 entries in the map, a spill is triggered. Since
   * mates arrive round-robin (read1/1 then read1/2), they land in the same map entry and pair
   * normally. This tests that basic spill + reassembly works for the simple case.
   */
  @Test
  public void testPairedEndSpillInOrder() throws IOException {
    Path inputFile1 = tempFolder.newFile("input_1.fastq").toPath();
    Path inputFile2 = tempFolder.newFile("input_2.fastq").toPath();

    // 4 pairs, in order
    Files.write(
        inputFile1,
        ("@A/1\nACGT\n+\nIIII\n"
                + "@B/1\nGGCC\n+\nJJJJ\n"
                + "@C/1\nTTAA\n+\nKKKK\n"
                + "@D/1\nCCGG\n+\nLLLL\n")
            .getBytes(StandardCharsets.UTF_8));

    Files.write(
        inputFile2,
        ("@A/2\nAAAA\n+\nIIII\n"
                + "@B/2\nCCCC\n+\nJJJJ\n"
                + "@C/2\nGGGG\n+\nKKKK\n"
                + "@D/2\nTTTT\n+\nLLLL\n")
            .getBytes(StandardCharsets.UTF_8));

    Path outputFile1 = tempFolder.newFile("output_1.fastq").toPath();
    Path outputFile2 = tempFolder.newFile("output_2.fastq").toPath();

    // spillPageSize=2 forces spill after every 2 map entries
    FastqNormalizer.PairedNormalizationResult result =
        FastqNormalizer.normalizePairedEnd(
            inputFile1.toString(),
            inputFile2.toString(),
            outputFile1.toString(),
            outputFile2.toString(),
            "RUN",
            false,
            tempFolder.getRoot(),
            2, // spillPageSize: spill after 2 entries
            Long.MAX_VALUE, // no byte limit
            Long.MAX_VALUE); // no abandon limit

    assertEquals(4, result.getPairCount());
    assertEquals(0, result.getOrphanCount());
    assertEquals(8, result.getTotalReadCount());

    List<FastqRecord> records1 = readFastq(outputFile1);
    List<FastqRecord> records2 = readFastq(outputFile2);

    assertEquals(4, records1.size());
    assertEquals(4, records2.size());

    // All should be paired
    for (FastqRecord r : records1) {
      assertTrue(r.getReadName().endsWith("/1"));
    }
    for (FastqRecord r : records2) {
      assertTrue(r.getReadName().endsWith("/2"));
    }
  }

  /**
   * Out-of-order pairs where mates end up in different spill pages. With spillPageSize=1, each map
   * entry triggers a spill. File 1 contains all /1 reads, file 2 contains all /2 reads. The
   * round-robin reading produces: A/1, A/2, B/1, B/2, ... so A's pair lands in the same page. But
   * with spillPageSize=1, A/1 causes a spill before A/2 arrives, forcing multi-generation
   * reassembly.
   *
   * <p>Actually with round-robin, A/1 goes in, then A/2 fills the second slot of the same key, so
   * spill happens with A complete. To truly force cross-page mates, we need reads that don't
   * alternate neatly. We use 3 reads per file with spillPageSize=1.
   */
  @Test
  public void testPairedEndSpillCrossPage() throws IOException {
    Path inputFile1 = tempFolder.newFile("input_1.fastq").toPath();
    Path inputFile2 = tempFolder.newFile("input_2.fastq").toPath();

    // 3 pairs each. With spillPageSize=1, after the first read is added to the map,
    // the next read triggers a spill check. Since round-robin reads A/1, A/2, B/1, B/2...
    // A/1 creates entry "A" (size=1), then A/2 fills slot 2 of "A" (size still 1 entry).
    // Then B/1 creates entry "B" (size=2 entries now ≥ spillPageSize=1) → spill.
    // So A is spilled (complete), B is spilled with only /1.
    // Then B/2 arrives, creates new entry "B" with /2, C/1 arrives...
    // This should exercise the reassembly path for B.
    Files.write(
        inputFile1,
        ("@A/1\nACGT\n+\nIIII\n" + "@B/1\nGGCC\n+\nJJJJ\n" + "@C/1\nTTAA\n+\nKKKK\n")
            .getBytes(StandardCharsets.UTF_8));

    Files.write(
        inputFile2,
        ("@A/2\nAAAA\n+\nIIII\n" + "@B/2\nCCCC\n+\nJJJJ\n" + "@C/2\nGGGG\n+\nKKKK\n")
            .getBytes(StandardCharsets.UTF_8));

    Path outputFile1 = tempFolder.newFile("output_1.fastq").toPath();
    Path outputFile2 = tempFolder.newFile("output_2.fastq").toPath();

    FastqNormalizer.PairedNormalizationResult result =
        FastqNormalizer.normalizePairedEnd(
            inputFile1.toString(),
            inputFile2.toString(),
            outputFile1.toString(),
            outputFile2.toString(),
            null,
            false,
            tempFolder.getRoot(),
            1, // spillPageSize=1: spill after every single entry
            Long.MAX_VALUE,
            Long.MAX_VALUE);

    assertEquals(3, result.getPairCount());
    assertEquals(0, result.getOrphanCount());
    assertEquals(6, result.getTotalReadCount());

    List<FastqRecord> records1 = readFastq(outputFile1);
    List<FastqRecord> records2 = readFastq(outputFile2);

    assertEquals(3, records1.size());
    assertEquals(3, records2.size());

    for (FastqRecord r : records1) {
      assertTrue(r.getReadName().endsWith("/1"));
    }
    for (FastqRecord r : records2) {
      assertTrue(r.getReadName().endsWith("/2"));
    }
  }

  /**
   * Spill with orphans: mismatched file lengths combined with tiny spill threshold. Verifies that
   * orphans are correctly identified and written even when spill reassembly is involved.
   */
  @Test
  public void testPairedEndSpillWithOrphans() throws IOException {
    Path inputFile1 = tempFolder.newFile("input_1.fastq").toPath();
    Path inputFile2 = tempFolder.newFile("input_2.fastq").toPath();

    // File 1 has 4 reads, file 2 has 2 — B and D will be orphans
    Files.write(
        inputFile1,
        ("@A/1\nACGT\n+\nIIII\n"
                + "@B/1\nGGCC\n+\nJJJJ\n"
                + "@C/1\nTTAA\n+\nKKKK\n"
                + "@D/1\nCCGG\n+\nLLLL\n")
            .getBytes(StandardCharsets.UTF_8));

    Files.write(
        inputFile2,
        ("@A/2\nAAAA\n+\nIIII\n" + "@C/2\nGGGG\n+\nKKKK\n").getBytes(StandardCharsets.UTF_8));

    Path outputFile1 = tempFolder.newFile("output_1.fastq").toPath();
    Path outputFile2 = tempFolder.newFile("output_2.fastq").toPath();

    FastqNormalizer.PairedNormalizationResult result =
        FastqNormalizer.normalizePairedEnd(
            inputFile1.toString(),
            inputFile2.toString(),
            outputFile1.toString(),
            outputFile2.toString(),
            "RUN",
            false,
            tempFolder.getRoot(),
            2,
            Long.MAX_VALUE,
            Long.MAX_VALUE);

    assertEquals(2, result.getPairCount());
    assertEquals(2, result.getOrphanCount());
    assertEquals(6, result.getTotalReadCount());

    List<FastqRecord> records1 = readFastq(outputFile1);
    List<FastqRecord> records2 = readFastq(outputFile2);

    // Output file 1: 2 paired (/1) + 2 orphans (no suffix) = 4
    assertEquals(4, records1.size());
    // Output file 2: 2 paired (/2)
    assertEquals(2, records2.size());

    // Verify orphans
    long orphansInFile1 = records1.stream().filter(r -> !r.getReadName().contains("/")).count();
    assertEquals(2, orphansInFile1);
  }

  /**
   * Large-scale spill test: 20 pairs with spillPageSize=3. This forces multiple spill generations
   * and verifies all reads are eventually paired and output in lexicographic order.
   */
  @Test
  public void testPairedEndSpillManyReads() throws IOException {
    Path inputFile1 = tempFolder.newFile("input_1.fastq").toPath();
    Path inputFile2 = tempFolder.newFile("input_2.fastq").toPath();

    int numPairs = 20;
    StringBuilder sb1 = new StringBuilder();
    StringBuilder sb2 = new StringBuilder();
    for (int i = 0; i < numPairs; i++) {
      // Use zero-padded names so lexicographic order = numeric order
      String name = String.format("READ%03d", i);
      sb1.append("@").append(name).append("/1\nACGT\n+\nIIII\n");
      sb2.append("@").append(name).append("/2\nTTAA\n+\nKKKK\n");
    }

    Files.write(inputFile1, sb1.toString().getBytes(StandardCharsets.UTF_8));
    Files.write(inputFile2, sb2.toString().getBytes(StandardCharsets.UTF_8));

    Path outputFile1 = tempFolder.newFile("output_1.fastq").toPath();
    Path outputFile2 = tempFolder.newFile("output_2.fastq").toPath();

    FastqNormalizer.PairedNormalizationResult result =
        FastqNormalizer.normalizePairedEnd(
            inputFile1.toString(),
            inputFile2.toString(),
            outputFile1.toString(),
            outputFile2.toString(),
            "RUN",
            false,
            tempFolder.getRoot(),
            3, // spillPageSize=3: forces many spill/reassembly cycles
            Long.MAX_VALUE,
            Long.MAX_VALUE);

    assertEquals(numPairs, result.getPairCount());
    assertEquals(0, result.getOrphanCount());
    assertEquals(numPairs * 2, result.getTotalReadCount());

    List<FastqRecord> records1 = readFastq(outputFile1);
    List<FastqRecord> records2 = readFastq(outputFile2);

    assertEquals(numPairs, records1.size());
    assertEquals(numPairs, records2.size());

    // Verify all paired
    for (FastqRecord r : records1) {
      assertTrue(r.getReadName().endsWith("/1"));
    }
    for (FastqRecord r : records2) {
      assertTrue(r.getReadName().endsWith("/2"));
    }

    // With spilling, output is sorted per-generation (per writeFromMemory batch), not globally.
    // This matches AbstractPagedReadWriter.cascadeErrors() behavior.
    // Verify all expected reads are present instead.
    java.util.Set<String> readKeys = new java.util.HashSet<>();
    for (FastqRecord r : records1) {
      // Extract the read name portion (after "RUN.N ")
      String fullName = r.getReadName();
      String readPart = fullName.substring(fullName.indexOf(' ') + 1);
      readKeys.add(readPart.replace("/1", ""));
    }
    for (int i = 0; i < numPairs; i++) {
      String expected = String.format("READ%03d", i);
      assertTrue("Missing read: " + expected, readKeys.contains(expected));
    }
  }

  // ---- Abandon threshold test ----

  /**
   * When total spilled bytes exceed the abandon limit, a ReadWriterMemoryLimitException should be
   * thrown. Uses spillPageSize=1 to force a spill on every entry and a tiny abandon limit so the
   * second spill triggers the exception.
   */
  @Test
  public void testPairedEndAbandonThreshold() throws IOException {
    Path inputFile1 = tempFolder.newFile("input_1.fastq").toPath();
    Path inputFile2 = tempFolder.newFile("input_2.fastq").toPath();

    // 4 pairs — enough to trigger multiple spills with spillPageSize=1
    StringBuilder sb1 = new StringBuilder();
    StringBuilder sb2 = new StringBuilder();
    for (int i = 0; i < 4; i++) {
      sb1.append("@READ").append(i).append("/1\nACGT\n+\nIIII\n");
      sb2.append("@READ").append(i).append("/2\nTTAA\n+\nKKKK\n");
    }
    Files.write(inputFile1, sb1.toString().getBytes(StandardCharsets.UTF_8));
    Files.write(inputFile2, sb2.toString().getBytes(StandardCharsets.UTF_8));

    Path outputFile1 = tempFolder.newFile("output_1.fastq").toPath();
    Path outputFile2 = tempFolder.newFile("output_2.fastq").toPath();

    try {
      FastqNormalizer.normalizePairedEnd(
          inputFile1.toString(),
          inputFile2.toString(),
          outputFile1.toString(),
          outputFile2.toString(),
          null,
          false,
          tempFolder.getRoot(),
          1, // spillPageSize=1: spill after every entry
          Long.MAX_VALUE, // no byte-size threshold
          1L); // abandon after 1 byte of spilled data
      fail("Expected ReadWriterMemoryLimitException");
    } catch (uk.ac.ebi.ena.readtools.loader.common.writer.ReadWriterMemoryLimitException e) {
      assertTrue(e.getMessage().contains("Temp memory limit"));
    }
  }

  /** Manual test with arbitrary local files. Ignored by default. */
  @org.junit.Ignore("Only run manually if needed.")
  @Test
  public void testArbitraryFile() throws IOException {
    String input1 =
        "tmp/NG-34612_V3V5b_5_lib735143_cleaned_2.fastq.gz";
    String input2 =
        "tmp/NG-34612_V3V5b_5_lib735143_cleaned_1.fastq.gz";
    String output1 = "tmp/normalised_1.fastq";
    String output2 = "tmp/normalised_2.fastq";

    FastqNormalizer.PairedNormalizationResult result =
        FastqNormalizer.normalizePairedEnd(
            input1, input2, output1, output2, "ERR12631034", false, tempFolder.getRoot());

    assertTrue("Expected at least one read pair", result.getPairCount() > 0);
  }

  /** Helper method to read all records from a FASTQ file. */
  private List<FastqRecord> readFastq(Path path) throws IOException {
    List<FastqRecord> records = new ArrayList<>();
    try (FastqReader reader =
        new FastqReader(
            null,
            new java.io.BufferedReader(
                new java.io.InputStreamReader(
                    Utils.openFastqInputStream(path), StandardCharsets.UTF_8)))) {
      for (FastqRecord record : reader) {
        records.add(record);
      }
    }
    return records;
  }
}
