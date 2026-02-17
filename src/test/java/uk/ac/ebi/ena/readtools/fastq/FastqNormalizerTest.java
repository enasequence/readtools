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
import org.junit.Ignore;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;
import uk.ac.ebi.ena.readtools.utils.Utils;

public class FastqNormalizerTest {

  @Rule public TemporaryFolder tempFolder = new TemporaryFolder();

  @Test
  public void testSingleEndSangerPassthrough() throws IOException {
    // Create a simple Sanger FASTQ file
    // Use quality scores that are unambiguous Sanger (ASCII < 64)
    Path inputFile = tempFolder.newFile("input.fastq").toPath();
    Files.write(
        inputFile,
        ("@READ1\n"
                + "ACGT\n"
                + "+\n"
                + "###$\n" // Phred+33 quality scores (ASCII 35,35,35,36 = Phred 2,2,2,3)
                + "@READ2\n"
                + "GGCC\n"
                + "+\n"
                + "%%&&\n") // Phred+33 quality scores (ASCII 37,37,38,38 = Phred 4,4,5,5)
            .getBytes(StandardCharsets.UTF_8));

    Path outputFile = tempFolder.newFile("output.fastq").toPath();

    // Normalize without prefix
    long count =
        FastqNormalizer.normalizeSingleEnd(
            inputFile.toString(), outputFile.toString(), null, false);

    assertEquals(2, count);

    // Read and verify output
    List<FastqRecord> records = readFastq(outputFile);
    assertEquals(2, records.size());

    // Verify read 1
    assertEquals("READ1", records.get(0).getReadName());
    assertEquals("ACGT", records.get(0).getReadString());
    assertEquals("###$", records.get(0).getBaseQualityString());

    // Verify read 2
    assertEquals("READ2", records.get(1).getReadName());
    assertEquals("GGCC", records.get(1).getReadString());
    assertEquals("%%&&", records.get(1).getBaseQualityString());
  }

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

    // Normalize with prefix
    long count =
        FastqNormalizer.normalizeSingleEnd(
            inputFile.toString(), outputFile.toString(), "ERR123456", false);

    assertEquals(3, count);

    // Read and verify output
    List<FastqRecord> records = readFastq(outputFile);
    assertEquals(3, records.size());

    // Verify prefix format: {prefix}.{counter} {originalName}
    assertEquals("ERR123456.1 READ1", records.get(0).getReadName());
    assertEquals("ERR123456.2 READ2", records.get(1).getReadName());
    assertEquals("ERR123456.3 READ3", records.get(2).getReadName());
  }

  @Test
  public void testSingleEndWithoutPrefix() throws IOException {
    Path inputFile = tempFolder.newFile("input.fastq").toPath();
    Files.write(
        inputFile, ("@READ1\n" + "ACGT\n" + "+\n" + "IIII\n").getBytes(StandardCharsets.UTF_8));

    Path outputFile = tempFolder.newFile("output.fastq").toPath();

    // Normalize without prefix
    long count =
        FastqNormalizer.normalizeSingleEnd(
            inputFile.toString(), outputFile.toString(), null, false);

    assertEquals(1, count);

    // Debug: check if file exists and has content
    assertTrue("Output file should exist", Files.exists(outputFile));
    long fileSize = Files.size(outputFile);
    assertTrue("Output file should not be empty, size=" + fileSize, fileSize > 0);

    // Debug: print first few lines
    System.out.println("=== Output file content ===");
    Files.lines(outputFile).limit(10).forEach(System.out::println);
    System.out.println("=== End output ===");

    // Read and verify output
    List<FastqRecord> records = readFastq(outputFile);
    assertEquals(1, records.size());
    assertEquals("READ1", records.get(0).getReadName());
  }

  @Test
  public void testSingleEndUracilConversion() throws IOException {
    Path inputFile = tempFolder.newFile("input.fastq").toPath();
    Files.write(
        inputFile,
        ("@READ1\n"
                + "AUGCU\n" // Mix of U and T
                + "+\n"
                + "IIIII\n"
                + "@READ2\n"
                + "gacggaUCuauagcaaaacu\n" // Mix of upper and lower case U
                + "+\n"
                + "IIIIIIIIIIIIIIIIIIII\n")
            .getBytes(StandardCharsets.UTF_8));

    Path outputFile = tempFolder.newFile("output.fastq").toPath();

    // Normalize with uracil conversion
    long count =
        FastqNormalizer.normalizeSingleEnd(inputFile.toString(), outputFile.toString(), null, true);

    assertEquals(2, count);

    // Read and verify output
    List<FastqRecord> records = readFastq(outputFile);
    assertEquals(2, records.size());

    // Verify U→T conversion (lowercase u→t, uppercase U→T)
    assertEquals("ATGCT", records.get(0).getReadString());
    assertEquals("gacggaTCTaTagcaaaacT", records.get(1).getReadString());
  }

  @Test
  public void testPairedEndInOrder() throws IOException {
    // Create paired FASTQ files with reads in matching order
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

    // Normalize paired-end
    long count =
        FastqNormalizer.normalizePairedEnd(
            inputFile1.toString(),
            inputFile2.toString(),
            outputFile1.toString(),
            outputFile2.toString(),
            null,
            false,
            tempFolder.getRoot());

    assertEquals(2, count);

    // Read and verify both output files
    List<FastqRecord> records1 = readFastq(outputFile1);
    List<FastqRecord> records2 = readFastq(outputFile2);

    assertEquals(2, records1.size());
    assertEquals(2, records2.size());

    // Verify pairing and /1 /2 suffixes
    assertTrue(records1.get(0).getReadName().endsWith("/1"));
    assertTrue(records2.get(0).getReadName().endsWith("/2"));
    assertTrue(records1.get(1).getReadName().endsWith("/1"));
    assertTrue(records2.get(1).getReadName().endsWith("/2"));
  }

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

    // Normalize with prefix
    long count =
        FastqNormalizer.normalizePairedEnd(
            inputFile1.toString(),
            inputFile2.toString(),
            outputFile1.toString(),
            outputFile2.toString(),
            "SRR654321",
            false,
            tempFolder.getRoot());

    assertEquals(1, count);

    // Read and verify output
    List<FastqRecord> records1 = readFastq(outputFile1);
    List<FastqRecord> records2 = readFastq(outputFile2);

    // Verify prefix format: {prefix}.{counter} {originalName}/1
    assertEquals("SRR654321.1 READ1/1", records1.get(0).getReadName());
    assertEquals("SRR654321.1 READ1/2", records2.get(0).getReadName());
  }

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

    // Mismatched lengths are tolerated: READ1 pairs normally, READ2 becomes an orphan
    long count =
        FastqNormalizer.normalizePairedEnd(
            inputFile1.toString(),
            inputFile2.toString(),
            outputFile1.toString(),
            outputFile2.toString(),
            null,
            false,
            tempFolder.getRoot());

    // 1 complete pair + 1 orphan = 2 total
    assertEquals(2, count);

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

  @Test
  public void testPairedEndCasavaFormat() throws IOException {
    // Test with Casava 1.8 format reads
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

    // Normalize
    long count =
        FastqNormalizer.normalizePairedEnd(
            inputFile1.toString(),
            inputFile2.toString(),
            outputFile1.toString(),
            outputFile2.toString(),
            null,
            false,
            tempFolder.getRoot());

    assertEquals(1, count);

    // Verify output exists and is valid
    List<FastqRecord> records1 = readFastq(outputFile1);
    List<FastqRecord> records2 = readFastq(outputFile2);

    assertEquals(1, records1.size());
    assertEquals(1, records2.size());
  }

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

    // Normalize with uracil conversion
    long count =
        FastqNormalizer.normalizePairedEnd(
            inputFile1.toString(),
            inputFile2.toString(),
            outputFile1.toString(),
            outputFile2.toString(),
            null,
            true,
            tempFolder.getRoot());

    assertEquals(1, count);

    // Read and verify output
    List<FastqRecord> records1 = readFastq(outputFile1);
    List<FastqRecord> records2 = readFastq(outputFile2);

    // Verify U→T conversion
    assertEquals("ATGCT", records1.get(0).getReadString());
    assertEquals("ttAAtt", records2.get(0).getReadString());
  }

  @Test
  public void testSingleEndRealFile() throws IOException {
    // Test with actual uracil-bases test file
    String inputFile = "src/test/resources/uracil-bases_1.fastq";
    Path outputFile = tempFolder.newFile("output.fastq").toPath();

    // Normalize with uracil conversion and prefix
    long count = FastqNormalizer.normalizeSingleEnd(inputFile, outputFile.toString(), "TEST", true);

    assertTrue(count > 0);

    // Read and verify output
    List<FastqRecord> records = readFastq(outputFile);
    assertEquals(count, records.size());

    // Verify prefix is applied
    assertTrue(records.get(0).getReadName().startsWith("TEST.1 "));

    // Verify no U bases remain
    for (FastqRecord record : records) {
      assertFalse(record.getReadString().contains("U"));
      assertFalse(record.getReadString().contains("u"));
    }
  }

  @Test
  public void testPairedEndRealFiles() throws IOException {
    // Test with actual paired FASTQ files
    String inputFile1 = "src/test/resources/2fastq/28239_1822_1.fastq";
    String inputFile2 = "src/test/resources/2fastq/28239_1822_2.fastq";
    Path outputFile1 = tempFolder.newFile("output_1.fastq").toPath();
    Path outputFile2 = tempFolder.newFile("output_2.fastq").toPath();

    // Normalize
    long count =
        FastqNormalizer.normalizePairedEnd(
            inputFile1,
            inputFile2,
            outputFile1.toString(),
            outputFile2.toString(),
            "RUN001",
            false,
            tempFolder.getRoot());

    assertTrue(count > 0);

    // Read and verify output
    List<FastqRecord> records1 = readFastq(outputFile1);
    List<FastqRecord> records2 = readFastq(outputFile2);

    assertEquals(count, records1.size());
    assertEquals(count, records2.size());

    // Verify prefix and /1 /2 suffixes
    assertTrue(records1.get(0).getReadName().startsWith("RUN001.1 "));
    assertTrue(records1.get(0).getReadName().endsWith("/1"));
    assertTrue(records2.get(0).getReadName().startsWith("RUN001.1 "));
    assertTrue(records2.get(0).getReadName().endsWith("/2"));
  }

  @Ignore("Only run manually if needed.")
  @Test
  public void testArbitraryFile() throws IOException {
    String input1 = "/path/input_2.fastq.gz";
    String input2 = "/path/input_1.fastq.gz";
    String output1 = "/path/normalised_1.fastq";
    String output2 = "/path/normalised_2.fastq";

    long count =
        FastqNormalizer.normalizePairedEnd(
            input1, input2, output1, output2, "ERR12631034", false, tempFolder.getRoot());

    assertTrue("Expected at least one read pair", count > 0);
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
