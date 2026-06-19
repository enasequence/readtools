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
package uk.ac.ebi.ena.readtools.fastq.ena;

import static junit.framework.Assert.assertTrue;
import static junit.framework.TestCase.fail;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import org.apache.commons.lang3.exception.ExceptionUtils;
import org.junit.Assert;
import org.junit.Ignore;
import org.junit.Test;
import uk.ac.ebi.ena.readtools.common.reads.QualityNormalizer;
import uk.ac.ebi.ena.readtools.loader.common.FileCompression;
import uk.ac.ebi.ena.readtools.loader.common.InvalidBaseCharacterException;
import uk.ac.ebi.ena.readtools.loader.common.converter.ConverterException;
import uk.ac.ebi.ena.readtools.loader.common.writer.ReadWriterException;
import uk.ac.ebi.ena.readtools.loader.common.writer.ReadWriterMemoryLimitException;
import uk.ac.ebi.ena.readtools.sam.Sam2Fastq;
import uk.ac.ebi.ena.readtools.utils.Utils;

/**
 * TODO The tests in here need to compare data between source fastq(s) and generated bam file as
 * well.
 */
public class Fastq2SamTest {
  @Test
  public void singleFastqReadAndBaseCount()
      throws IOException, ConverterException, ReadWriterException {

    Fastq2Sam.Params params = new Fastq2Sam.Params();
    params.tmp_root = System.getProperty("java.io.tmpdir");
    params.sample_name = "SM-001";
    params.data_file = Files.createTempFile(null, ".bam").toString();
    params.compression = FileCompression.NONE.name();
    params.files =
        Arrays.asList(
            new File(
                    Fastq2SamTest.class
                        .getClassLoader()
                        .getResource("fastq_spots_correct_paired_with_unpaired_1.txt")
                        .getFile())
                .getAbsolutePath());

    Fastq2Sam fastq2Sam = new Fastq2Sam();
    fastq2Sam.create(params);

    Assert.assertTrue(new File(params.data_file).length() > 0);
    Assert.assertEquals(4, fastq2Sam.getTotalReadCount());
    Assert.assertEquals(404, fastq2Sam.getTotalBaseCount());
    assertBamMatchesFastqInputs(params);
  }

  @Test
  public void pairedFastqReadAndBaseCount()
      throws IOException, ConverterException, ReadWriterException {

    Fastq2Sam.Params params = new Fastq2Sam.Params();
    params.tmp_root = System.getProperty("java.io.tmpdir");
    params.sample_name = "SM-001";
    params.data_file = Files.createTempFile(null, ".bam").toString();
    params.compression = FileCompression.NONE.name();
    params.files =
        Arrays.asList(
            new File(
                    Fastq2SamTest.class
                        .getClassLoader()
                        .getResource("fastq_spots_correct_paired_with_unpaired_1.txt")
                        .getFile())
                .getAbsolutePath(),
            new File(
                    Fastq2SamTest.class
                        .getClassLoader()
                        .getResource("fastq_spots_correct_paired_with_unpaired_2.txt")
                        .getFile())
                .getAbsolutePath());

    Fastq2Sam fastq2Sam = new Fastq2Sam();
    fastq2Sam.create(params);

    Assert.assertTrue(new File(params.data_file).length() > 0);
    Assert.assertEquals(8, fastq2Sam.getTotalReadCount());
    Assert.assertEquals(808, fastq2Sam.getTotalBaseCount());
    assertBamMatchesFastqInputs(params);
  }

  @Test
  public void pairedFastqCasavaLikeNoPairNum()
      throws IOException, ConverterException, ReadWriterException {

    Fastq2Sam.Params params = new Fastq2Sam.Params();
    params.tmp_root = System.getProperty("java.io.tmpdir");
    params.sample_name = "SM-001";
    params.data_file = Files.createTempFile(null, ".bam").toString();
    params.compression = FileCompression.NONE.name();
    params.files =
        Arrays.asList(
            new File(
                    Fastq2SamTest.class
                        .getClassLoader()
                        .getResource("fastq_spots_correct2_1b.txt")
                        .getFile())
                .getAbsolutePath(),
            new File(
                    Fastq2SamTest.class
                        .getClassLoader()
                        .getResource("fastq_spots_correct2_2b.txt")
                        .getFile())
                .getAbsolutePath());

    Fastq2Sam fastq2Sam = new Fastq2Sam();
    fastq2Sam.create(params);

    Assert.assertTrue(new File(params.data_file).length() > 0);
    Assert.assertEquals(6, fastq2Sam.getTotalReadCount());
    Assert.assertEquals(906, fastq2Sam.getTotalBaseCount());
    assertBamMatchesFastqInputs(params);
  }

  @Test
  public void pairedFastqCasavaLikeNoPairNumFlowcellDash()
      throws IOException, ConverterException, ReadWriterException {

    Fastq2Sam.Params params = new Fastq2Sam.Params();
    params.tmp_root = System.getProperty("java.io.tmpdir");
    params.sample_name = "SM-001";
    params.data_file = Files.createTempFile(null, ".bam").toString();
    params.compression = FileCompression.NONE.name();
    params.files =
        Arrays.asList(
            new File(
                    Fastq2SamTest.class
                        .getClassLoader()
                        .getResource("fastq_spots_correct2_1d.txt")
                        .getFile())
                .getAbsolutePath(),
            new File(
                    Fastq2SamTest.class
                        .getClassLoader()
                        .getResource("fastq_spots_correct2_2d.txt")
                        .getFile())
                .getAbsolutePath());

    Fastq2Sam fastq2Sam = new Fastq2Sam();
    fastq2Sam.create(params);

    Assert.assertTrue(new File(params.data_file).length() > 0);
    Assert.assertEquals(6, fastq2Sam.getTotalReadCount());
    Assert.assertEquals(906, fastq2Sam.getTotalBaseCount());
    assertBamMatchesFastqInputs(params);
  }

  @Test
  public void pairedFastqCasavaLikeSingleDigitYNoPairNum()
      throws IOException, ConverterException, ReadWriterException {

    Fastq2Sam.Params params = new Fastq2Sam.Params();
    params.tmp_root = System.getProperty("java.io.tmpdir");
    params.sample_name = "SM-001";
    params.data_file = Files.createTempFile(null, ".bam").toString();
    params.compression = FileCompression.NONE.name();
    params.files =
        Arrays.asList(
            new File(
                    Fastq2SamTest.class
                        .getClassLoader()
                        .getResource("fastq_spots_correct2_1c.txt")
                        .getFile())
                .getAbsolutePath(),
            new File(
                    Fastq2SamTest.class
                        .getClassLoader()
                        .getResource("fastq_spots_correct2_2c.txt")
                        .getFile())
                .getAbsolutePath());

    Fastq2Sam fastq2Sam = new Fastq2Sam();
    fastq2Sam.create(params);

    Assert.assertTrue(new File(params.data_file).length() > 0);
    Assert.assertEquals(6, fastq2Sam.getTotalReadCount());
    Assert.assertEquals(906, fastq2Sam.getTotalBaseCount());
    assertBamMatchesFastqInputs(params);
  }

  @Ignore("Only run manually if needed.")
  @Test
  public void twoLargeFilesManualTest()
      throws IOException, ConverterException, ReadWriterException {
    Fastq2Sam.Params params = new Fastq2Sam.Params();
    params.tmp_root = System.getProperty("java.io.tmpdir");
    params.sample_name = "SM-001";
    params.data_file = Files.createTempFile(null, ".bam").toString();
    params.compression = FileCompression.GZ.name();
    params.files =
        Arrays.asList(
            new File(Fastq2SamTest.class.getClassLoader().getResource("F1.fastq.gz").getFile())
                .getAbsolutePath(),
            new File(Fastq2SamTest.class.getClassLoader().getResource("F2.fastq.gz").getFile())
                .getAbsolutePath());

    params.spill_page_size_bytes = 2L * 1024L * 1024L * 1024L;
    params.spill_abandon_limit_bytes = 10L * 1024L * 1024L * 1024L;

    Fastq2Sam fastq2Sam = new Fastq2Sam();
    fastq2Sam.create(params);

    Assert.assertTrue(new File(params.data_file).length() > 0);
  }

  @Test
  public void twoReadPairNumberDetection()
      throws IOException, ConverterException, ReadWriterException {
    Fastq2Sam.Params params = new Fastq2Sam.Params();
    params.tmp_root = System.getProperty("java.io.tmpdir");
    params.sample_name = "SM-001";
    params.data_file = Files.createTempFile(null, ".bam").toString();
    params.compression = FileCompression.GZ.name();
    params.files =
        Arrays.asList(
            new File(
                    Fastq2SamTest.class
                        .getClassLoader()
                        .getResource("read_pair_number_detection_1.fastq.gz")
                        .getFile())
                .getAbsolutePath(),
            new File(
                    Fastq2SamTest.class
                        .getClassLoader()
                        .getResource("read_pair_number_detection_2.fastq.gz")
                        .getFile())
                .getAbsolutePath());

    params.spill_page_size_bytes = 2L * 1024L * 1024L * 1024L;
    params.spill_abandon_limit_bytes = 10L * 1024L * 1024L * 1024L;

    Fastq2Sam fastq2Sam = new Fastq2Sam();
    fastq2Sam.create(params);

    Assert.assertTrue(new File(params.data_file).length() > 0);
  }

  @Test
  public void pairedFastqPairNumber3() throws IOException, ConverterException, ReadWriterException {

    Fastq2Sam.Params params = new Fastq2Sam.Params();
    params.tmp_root = System.getProperty("java.io.tmpdir");
    params.sample_name = "SM-001";
    params.data_file = Files.createTempFile(null, ".bam").toString();
    params.compression = FileCompression.NONE.name();
    params.files =
        Arrays.asList(
            new File(
                    Fastq2SamTest.class
                        .getClassLoader()
                        .getResource("fastq_spots_read_num_eq3_1.fastq")
                        .getFile())
                .getAbsolutePath(),
            new File(
                    Fastq2SamTest.class
                        .getClassLoader()
                        .getResource("fastq_spots_read_num_eq3_2.fastq")
                        .getFile())
                .getAbsolutePath());

    Fastq2Sam fastq2Sam = new Fastq2Sam();
    fastq2Sam.create(params);

    Assert.assertTrue(new File(params.data_file).length() > 0);
    Assert.assertEquals(6, fastq2Sam.getTotalReadCount());
    Assert.assertEquals(906, fastq2Sam.getTotalBaseCount());
    assertBamMatchesFastqInputs(params);
  }

  @Test
  public void testCasava18BarcodePreservedInBam()
      throws IOException, ConverterException, ReadWriterException {
    Path tempDir = Files.createTempDirectory("fastq_casava_test");

    Path fastqFile1 = tempDir.resolve("casava_1.fastq");
    Path fastqFile2 = tempDir.resolve("casava_2.fastq");

    // Read 1: filter=Y (failed), barcode=ATCACG+GCGCTA
    // Read 2: filter=N (passed), barcode=TTAGGC
    Files.write(
        fastqFile1,
        Arrays.asList(
            "@EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG+GCGCTA",
            "GATCGGAAGAGCACACGTCTGAACTCCAGTCA",
            "+",
            "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF",
            "@EAS139:136:FC706VJ:2:2104:15343:197394 1:N:0:TTAGGC",
            "GATCGGAAGAGCACACGTCTGAACTCCAGTCA",
            "+",
            "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"));

    Files.write(
        fastqFile2,
        Arrays.asList(
            "@EAS139:136:FC706VJ:2:2104:15343:197393 2:Y:18:ATCACG+GCGCTA",
            "GATCGGAAGAGCACACGTCTGAACTCCAGTCA",
            "+",
            "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF",
            "@EAS139:136:FC706VJ:2:2104:15343:197394 2:N:0:TTAGGC",
            "GATCGGAAGAGCACACGTCTGAACTCCAGTCA",
            "+",
            "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"));

    Fastq2Sam.Params params = new Fastq2Sam.Params();
    params.tmp_root = tempDir.toString();
    params.sample_name = "SM-001";
    params.data_file = Files.createTempFile(tempDir, "output", ".bam").toString();
    params.compression = FileCompression.NONE.name();
    params.files = Arrays.asList(fastqFile1.toString(), fastqFile2.toString());

    Fastq2Sam fastq2Sam = new Fastq2Sam();
    fastq2Sam.create(params);

    Assert.assertEquals(4, fastq2Sam.getTotalReadCount());

    try (SamReader samReader = SamReaderFactory.makeDefault().open(new File(params.data_file))) {
      for (final SAMRecord rec : samReader) {
        String qname = rec.getReadName();
        String bc = (String) rec.getAttribute("BC");

        if (qname.endsWith("197393")) {
          Assert.assertEquals("ATCACG+GCGCTA", bc);
          Assert.assertTrue(
              "QCFAIL should be set for filter=Y", rec.getReadFailsVendorQualityCheckFlag());
        } else if (qname.endsWith("197394")) {
          Assert.assertEquals("TTAGGC", bc);
          Assert.assertFalse(
              "QCFAIL should not be set for filter=N", rec.getReadFailsVendorQualityCheckFlag());
        } else {
          Assert.fail("Unexpected QNAME: " + qname);
        }
      }
    }
  }

  @Test
  public void testMemoryPaging() throws IOException, ConverterException, ReadWriterException {

    Fastq2Sam.Params params = new Fastq2Sam.Params();
    params.tmp_root = System.getProperty("java.io.tmpdir");
    params.sample_name = "SM-001";
    params.data_file = Files.createTempFile(null, ".bam").toString();
    params.compression = FileCompression.NONE.name();
    params.files =
        Arrays.asList(
            new File(
                    Fastq2SamTest.class
                        .getClassLoader()
                        .getResource("fastq_spots_correct_paired_with_unpaired_1.txt")
                        .getFile())
                .getAbsolutePath(),
            new File(
                    Fastq2SamTest.class
                        .getClassLoader()
                        .getResource("fastq_spots_correct_paired_with_unpaired_2.txt")
                        .getFile())
                .getAbsolutePath());

    params.spill_page_size_bytes = 400L;
    params.spill_abandon_limit_bytes = 10L * 1024L;

    Fastq2Sam fastq2Sam = new Fastq2Sam();
    fastq2Sam.create(params);

    Assert.assertTrue(new File(params.data_file).length() > 0);
    Assert.assertEquals(8, fastq2Sam.getTotalReadCount());
    Assert.assertEquals(808, fastq2Sam.getTotalBaseCount());
    assertBamMatchesFastqInputs(params);
  }

  @Test
  public void testMemoryPagingAbandonLimit()
      throws IOException, ConverterException, ReadWriterException {
    Fastq2Sam.Params params = new Fastq2Sam.Params();
    params.tmp_root = System.getProperty("java.io.tmpdir");
    params.sample_name = "SM-001";
    params.data_file = Files.createTempFile(null, ".bam").toString();
    params.compression = FileCompression.NONE.name();
    params.files =
        Arrays.asList(
            new File(
                    Fastq2SamTest.class
                        .getClassLoader()
                        .getResource("fastq_spots_correct_paired_with_unpaired_1a.txt")
                        .getFile())
                .getAbsolutePath(),
            new File(
                    Fastq2SamTest.class
                        .getClassLoader()
                        .getResource("fastq_spots_correct_paired_with_unpaired_2a.txt")
                        .getFile())
                .getAbsolutePath());

    params.spill_page_size_bytes = 16L;
    params.spill_abandon_limit_bytes = 16L;

    Fastq2Sam fastq2Sam = new Fastq2Sam();

    try {
      fastq2Sam.create(params);
      fail();
    } catch (ReadWriterMemoryLimitException e) {
      assertTrue(e.getMessage().contains("Temp memory limit"));
    }
  }

  @Test
  public void testCorrectPairedWithUnpaired()
      throws IOException, ConverterException, ReadWriterException {
    File inpFile1 =
        new File(
            Fastq2SamTest.class
                .getClassLoader()
                .getResource("fastq_spots_correct_paired_with_unpaired_1.txt")
                .getFile());
    File inpFile2 =
        new File(
            Fastq2SamTest.class
                .getClassLoader()
                .getResource("fastq_spots_correct_paired_with_unpaired_2.txt")
                .getFile());

    File outFile = Files.createTempFile(null, ".bam").toFile();

    Fastq2Sam.Params params = new Fastq2Sam.Params();
    params.tmp_root = System.getProperty("java.io.tmpdir");
    params.sample_name = "SM-001";
    params.data_file = outFile.getAbsolutePath();
    params.compression = FileCompression.NONE.name();
    params.files = Arrays.asList(inpFile1.getAbsolutePath(), inpFile2.getAbsolutePath());

    Fastq2Sam fastq2Sam = new Fastq2Sam();
    fastq2Sam.create(params);

    Assert.assertTrue(outFile.length() > 0);
    Assert.assertEquals(8, fastq2Sam.getTotalReadCount());
    Assert.assertEquals(808, fastq2Sam.getTotalBaseCount());

    Map<String, List<FastqRecord>> fastqRecordMap = createFastqRecordMap(inpFile1, inpFile2);

    int recCount = 0;

    try (SamReader samReader = SamReaderFactory.makeDefault().open(outFile); ) {
      Assert.assertEquals(
          params.sample_name, samReader.getFileHeader().getReadGroup("A").getSample());

      for (final SAMRecord rec : samReader) {
        ++recCount;

        List<FastqRecord> fastqRecs = fastqRecordMap.get(rec.getReadName());

        if (!rec.getReadPairedFlag()) {
          Assert.assertEquals(1, fastqRecs.size());
          Assert.assertEquals(fastqRecs.get(0).getReadString(), rec.getReadString());
        } else if (rec.getFirstOfPairFlag()) {
          Assert.assertEquals(fastqRecs.get(0).getReadString(), rec.getReadString());
        } else if (rec.getSecondOfPairFlag()) {
          Assert.assertEquals(fastqRecs.get(1).getReadString(), rec.getReadString());
        } else {
          Assert.fail("Unexpected pairing configuration of the SAM record.");
        }
      }
    }

    int expectedRecCount =
        fastqRecordMap.values().stream()
            .map(list -> list.size())
            .reduce((listSize1, listSize2) -> listSize1 + listSize2)
            .get();

    Assert.assertEquals(expectedRecCount, recCount);
  }

  @Test
  public void testUracilFastqSingle() throws IOException {
    Fastq2Sam.Params params = new Fastq2Sam.Params();
    params.tmp_root = System.getProperty("java.io.tmpdir");
    params.sample_name = "SM-001";
    params.compression = FileCompression.NONE.name();
    params.files =
        Arrays.asList(
            new File(
                    Fastq2SamTest.class
                        .getClassLoader()
                        .getResource("uracil-bases_1.fastq")
                        .getFile())
                .getAbsolutePath());

    // with flag set to false.
    params.convertUracil = false;
    params.data_file = Files.createTempFile(null, ".bam").toString();

    InvalidBaseCharacterException e = null;
    try {
      new Fastq2Sam().create(params);
    } catch (Exception ex) {
      e = (InvalidBaseCharacterException) ExceptionUtils.getRootCause(ex);
    }
    Assert.assertTrue(
        Pattern.matches(
            "[uU]{1,2}",
            e.getInvalidCharacters().stream().map(String::valueOf).collect(Collectors.joining())));

    // with flag set to true.
    params.convertUracil = true;
    params.data_file = Files.createTempFile(null, ".bam").toString();

    // Since U or u is not an acceptable base character inside a BAM file, execution of following
    // line without
    // any error would mean that Uracil base conversion did take place and there is no need for its
    // verification.
    Fastq2Sam fastq2Sam = new Fastq2Sam();
    fastq2Sam.create(params);

    Assert.assertTrue(new File(params.data_file).length() > 0);
    Assert.assertEquals(2, fastq2Sam.getTotalReadCount());
    Assert.assertEquals(40, fastq2Sam.getTotalBaseCount());
    assertBamMatchesFastqInputs(params);
  }

  @Test
  public void testUracilFastqPaired() throws IOException {
    Fastq2Sam.Params params = new Fastq2Sam.Params();
    params.tmp_root = System.getProperty("java.io.tmpdir");
    params.sample_name = "SM-001";
    params.compression = FileCompression.NONE.name();
    params.files =
        Arrays.asList(
            new File(
                    Fastq2SamTest.class
                        .getClassLoader()
                        .getResource("uracil-bases_1.fastq")
                        .getFile())
                .getAbsolutePath(),
            new File(
                    Fastq2SamTest.class
                        .getClassLoader()
                        .getResource("uracil-bases_2.fastq")
                        .getFile())
                .getAbsolutePath());

    // with flag set to false.
    params.convertUracil = false;
    params.data_file = Files.createTempFile(null, ".bam").toString();

    InvalidBaseCharacterException e = null;
    try {
      new Fastq2Sam().create(params);
    } catch (Exception ex) {
      e = (InvalidBaseCharacterException) ExceptionUtils.getRootCause(ex);
    }
    Assert.assertTrue(
        Pattern.matches(
            "[uU]{1,2}",
            e.getInvalidCharacters().stream().map(String::valueOf).collect(Collectors.joining())));

    // with flag set to true.
    params.convertUracil = true;
    params.data_file = Files.createTempFile(null, ".bam").toString();

    // Since U or u is not an acceptable base character inside a BAM file, execution of following
    // line without
    // any error would mean that Uracil base conversion did take place and there is no need for its
    // verification.
    Fastq2Sam fastq2Sam = new Fastq2Sam();
    fastq2Sam.create(params);

    Assert.assertTrue(new File(params.data_file).length() > 0);
    Assert.assertEquals(4, fastq2Sam.getTotalReadCount());
    Assert.assertEquals(80, fastq2Sam.getTotalBaseCount());
    assertBamMatchesFastqInputs(params);
  }

  @Test
  public void testInputFileOrderSwap() throws IOException {
    Path f1 = Paths.get("src/test/resources/ERR12387716/tf_1.fastq.gz");
    Path f2 = Paths.get("src/test/resources/ERR12387716/tf_2.fastq.gz");

    Fastq2Sam.Params params = new Fastq2Sam.Params();
    params.tmp_root = System.getProperty("java.io.tmpdir");
    params.sample_name = "SM-001";
    params.data_file = Files.createTempFile(null, ".bam").toString();
    params.compression = FileCompression.GZ.name();
    params.files = Arrays.asList(f1.toString(), f2.toString());

    params.spill_page_size_bytes = 2L * 1024L * 1024L * 1024L;
    params.spill_abandon_limit_bytes = 10L * 1024L * 1024L * 1024L;

    Fastq2Sam fastq2Sam1 = new Fastq2Sam();
    fastq2Sam1.create(params);

    Assert.assertTrue(new File(params.data_file).length() > 0);
    Assert.assertEquals(4, fastq2Sam1.getTotalReadCount());
    Assert.assertEquals(1203, fastq2Sam1.getTotalBaseCount());
    assertBamMatchesFastqInputs(params);
    List<String> firstRunRecords = canonicalBamRecords(new File(params.data_file), params);

    params.files = Arrays.asList(f2.toString(), f1.toString());
    Fastq2Sam fastq2Sam2 = new Fastq2Sam();
    fastq2Sam2.create(params);

    Assert.assertTrue(new File(params.data_file).length() > 0);
    Assert.assertEquals(4, fastq2Sam2.getTotalReadCount());
    Assert.assertEquals(1203, fastq2Sam2.getTotalBaseCount());
    assertBamMatchesFastqInputs(params);
    Assert.assertEquals(
        "Swapping input file order should not change logical BAM records",
        firstRunRecords,
        canonicalBamRecords(new File(params.data_file), params));
  }

  @Test
  public void testInvalidReadPairNumber() throws IOException {
    Path tempDir = Files.createTempDirectory("fastq_test");

    Path fastqFile1 = tempDir.resolve("test_invalid_1.fastq");
    Path fastqFile2 = tempDir.resolve("test_invalid_2.fastq");

    Files.write(
        fastqFile1,
        Arrays.asList(
            "@read1/1",
            "GATCGGAAGAGCACACGTCTGAACTCCAGTCA",
            "+",
            "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF",
            "@read2/3", // Invalid pair number
            "GATCGGAAGAGCACACGTCTGAACTCCAGTCA",
            "+",
            "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"));

    Files.write(
        fastqFile2,
        Arrays.asList(
            "@read1/2",
            "GATCGGAAGAGCACACGTCTGAACTCCAGTCA",
            "+",
            "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF",
            "@read2/4", // Invalid pair number
            "GATCGGAAGAGCACACGTCTGAACTCCAGTCA",
            "+",
            "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"));

    Fastq2Sam.Params params = new Fastq2Sam.Params();
    params.tmp_root = tempDir.toString();
    params.sample_name = "SM-001";
    params.data_file = Files.createTempFile(tempDir, "output", ".bam").toString();
    params.compression = FileCompression.NONE.name();
    params.files = Arrays.asList(fastqFile1.toString(), fastqFile2.toString());

    Fastq2Sam fastq2Sam = new Fastq2Sam();
    try {
      fastq2Sam.create(params);
      fail();
    } catch (ReadWriterException e) {
      Assert.assertTrue(e.getMessage().contains("Unexpected read pair number"));
    }
  }

  @Test
  public void testDuplicateReadSpot() throws IOException {
    Path tempDir = Files.createTempDirectory("fastq_test");

    Path fastqFile1 = tempDir.resolve("test_duplicate_1.fastq");
    Path fastqFile2 = tempDir.resolve("test_duplicate_2.fastq");

    Files.write(
        fastqFile1,
        Arrays.asList(
            "@read1/1",
            "GATCGGAAGAGCACACGTCTGAACTCCAGTCA",
            "+",
            "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"));

    Files.write(
        fastqFile2,
        Arrays.asList(
            "@read1/1", // Duplicate spot
            "GATCGGAAGAGCACACGTCTGAACTCCAGTCA",
            "+",
            "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF",
            "@read1/2",
            "GATCGGAAGAGCACACGTCTGAACTCCAGTCA",
            "+",
            "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"));

    Fastq2Sam.Params params = new Fastq2Sam.Params();
    params.tmp_root = tempDir.toString();
    params.sample_name = "SM-001";
    params.data_file = Files.createTempFile(tempDir, "output", ".bam").toString();
    params.compression = FileCompression.NONE.name();
    params.files = Arrays.asList(fastqFile1.toString(), fastqFile2.toString());

    Fastq2Sam fastq2Sam = new Fastq2Sam();
    try {
      fastq2Sam.create(params);
      fail();
    } catch (ReadWriterException e) {
      Assert.assertTrue(e.getMessage().contains("Got same spot twice"));
    }
  }

  @Test
  public void testCheckerboardPatternReads() throws IOException {
    Path tempDir = Files.createTempDirectory("fastq_test");

    Path fastqFile1 = tempDir.resolve("test_1.fastq");
    Path fastqFile2 = tempDir.resolve("test_2.fastq");

    Files.write(
        fastqFile1,
        Arrays.asList(
            "@read1/1",
            "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
            "+",
            "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF",
            "@read2/2", // Duplicate spot
            "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT",
            "+",
            "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF",
            "@read3/1",
            "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
            "+",
            "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF",
            "@read4/2",
            "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT",
            "+",
            "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"));

    Files.write(
        fastqFile2,
        Arrays.asList(
            "@read1/2",
            "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT",
            "+",
            "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF",
            "@read2/1",
            "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
            "+",
            "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF",
            "@read3/2",
            "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT",
            "+",
            "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF",
            "@read4/1",
            "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
            "+",
            "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"));

    Fastq2Sam.Params params = new Fastq2Sam.Params();
    params.tmp_root = tempDir.toString();
    params.sample_name = "SM-001";
    params.data_file = Files.createTempFile(tempDir, "output", ".bam").toString();
    params.compression = FileCompression.NONE.name();
    params.files = Arrays.asList(fastqFile1.toString(), fastqFile2.toString());

    Fastq2Sam fastq2Sam1 = new Fastq2Sam();
    fastq2Sam1.create(params);

    Assert.assertTrue(new File(params.data_file).length() > 0);
    Assert.assertEquals(8, fastq2Sam1.getTotalReadCount());
    Assert.assertEquals(256, fastq2Sam1.getTotalBaseCount());
    assertBamMatchesFastqInputs(params);
    List<String> firstRunRecords = canonicalBamRecords(new File(params.data_file), params);

    params.files = Arrays.asList(fastqFile2.toString(), fastqFile1.toString());
    Fastq2Sam fastq2Sam2 = new Fastq2Sam();
    fastq2Sam2.create(params);

    Assert.assertTrue(new File(params.data_file).length() > 0);
    Assert.assertEquals(8, fastq2Sam2.getTotalReadCount());
    Assert.assertEquals(256, fastq2Sam2.getTotalBaseCount());
    assertBamMatchesFastqInputs(params);
    Assert.assertEquals(
        "Swapping checkerboard input files should not change logical BAM records",
        firstRunRecords,
        canonicalBamRecords(new File(params.data_file), params));
  }

  // ---- Spill reassembly tests ----
  // These tests generate paired FASTQ files and use tiny spill_page_size values to force
  // disk spilling. They verify that multi-generation reassembly produces identical BAM output
  // to the no-spill case.

  /**
   * In-order pairs with spill_page_size=2. Mates arrive round-robin and land in the same map entry,
   * so pairs are complete before spilling. Verifies basic spill+reassembly produces identical BAM
   * to no-spill run.
   */
  @Test
  public void testSpillReassemblyInOrder()
      throws IOException, ConverterException, ReadWriterException {
    Path tempDir = Files.createTempDirectory("fastq_spill_test");

    Path fastqFile1 = tempDir.resolve("input_1.fastq");
    Path fastqFile2 = tempDir.resolve("input_2.fastq");

    // 6 pairs, in order
    List<String> lines1 = new ArrayList<>();
    List<String> lines2 = new ArrayList<>();
    for (int i = 0; i < 6; i++) {
      String name = String.format("READ%03d", i);
      lines1.addAll(
          Arrays.asList(
              "@" + name + "/1",
              "ACGTACGTACGTACGTACGTACGTACGTACGT",
              "+",
              "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"));
      lines2.addAll(
          Arrays.asList(
              "@" + name + "/2",
              "TTAATTAATTAATTAATTAATTAATTAATTAA",
              "+",
              "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"));
    }
    Files.write(fastqFile1, lines1);
    Files.write(fastqFile2, lines2);

    // Run without spilling (default thresholds)
    String noSpillBam = Files.createTempFile(tempDir, "nospill", ".bam").toString();
    Fastq2Sam.Params noSpillParams = new Fastq2Sam.Params();
    noSpillParams.tmp_root = tempDir.toString();
    noSpillParams.sample_name = "SM-001";
    noSpillParams.data_file = noSpillBam;
    noSpillParams.compression = FileCompression.NONE.name();
    noSpillParams.files = Arrays.asList(fastqFile1.toString(), fastqFile2.toString());

    Fastq2Sam noSpillF2s = new Fastq2Sam();
    noSpillF2s.create(noSpillParams);

    Assert.assertEquals(12, noSpillF2s.getTotalReadCount());
    assertBamMatchesFastqInputs(noSpillParams);
    List<String> noSpillRecords = canonicalBamRecords(new File(noSpillBam), noSpillParams);

    // Run with spill_page_size=2
    String spillBam = Files.createTempFile(tempDir, "spill", ".bam").toString();
    Fastq2Sam.Params spillParams = new Fastq2Sam.Params();
    spillParams.tmp_root = tempDir.toString();
    spillParams.sample_name = "SM-001";
    spillParams.data_file = spillBam;
    spillParams.compression = FileCompression.NONE.name();
    spillParams.files = Arrays.asList(fastqFile1.toString(), fastqFile2.toString());
    spillParams.spill_page_size = 2;
    spillParams.spill_page_size_bytes = Long.MAX_VALUE;
    spillParams.spill_abandon_limit_bytes = Long.MAX_VALUE;

    Fastq2Sam spillF2s = new Fastq2Sam();
    spillF2s.create(spillParams);

    Assert.assertEquals(12, spillF2s.getTotalReadCount());
    assertBamMatchesFastqInputs(spillParams);

    Assert.assertEquals(
        "Spill BAM should contain the same logical records as no-spill BAM",
        noSpillRecords,
        canonicalBamRecords(new File(spillBam), spillParams));
  }

  /**
   * Cross-page mates: with spill_page_size=1, mates from the same pair can end up in different
   * spill pages, requiring multi-generation reassembly. Verifies all reads are still correctly
   * paired and the BAM is identical to the no-spill baseline.
   */
  @Test
  public void testSpillReassemblyCrossPage()
      throws IOException, ConverterException, ReadWriterException {
    Path tempDir = Files.createTempDirectory("fastq_spill_test");

    Path fastqFile1 = tempDir.resolve("input_1.fastq");
    Path fastqFile2 = tempDir.resolve("input_2.fastq");

    // 4 pairs
    List<String> lines1 = new ArrayList<>();
    List<String> lines2 = new ArrayList<>();
    for (int i = 0; i < 4; i++) {
      String name = String.format("READ%03d", i);
      lines1.addAll(
          Arrays.asList(
              "@" + name + "/1",
              "ACGTACGTACGTACGTACGTACGTACGTACGT",
              "+",
              "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"));
      lines2.addAll(
          Arrays.asList(
              "@" + name + "/2",
              "TTAATTAATTAATTAATTAATTAATTAATTAA",
              "+",
              "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"));
    }
    Files.write(fastqFile1, lines1);
    Files.write(fastqFile2, lines2);

    // No-spill baseline
    String noSpillBam = Files.createTempFile(tempDir, "nospill", ".bam").toString();
    Fastq2Sam.Params noSpillParams = new Fastq2Sam.Params();
    noSpillParams.tmp_root = tempDir.toString();
    noSpillParams.sample_name = "SM-001";
    noSpillParams.data_file = noSpillBam;
    noSpillParams.compression = FileCompression.NONE.name();
    noSpillParams.files = Arrays.asList(fastqFile1.toString(), fastqFile2.toString());

    Fastq2Sam noSpillF2s = new Fastq2Sam();
    noSpillF2s.create(noSpillParams);

    Assert.assertEquals(8, noSpillF2s.getTotalReadCount());
    assertBamMatchesFastqInputs(noSpillParams);
    List<String> noSpillRecords = canonicalBamRecords(new File(noSpillBam), noSpillParams);

    // Spill with page_size=1 — forces spill after every single entry
    String spillBam = Files.createTempFile(tempDir, "spill", ".bam").toString();
    Fastq2Sam.Params spillParams = new Fastq2Sam.Params();
    spillParams.tmp_root = tempDir.toString();
    spillParams.sample_name = "SM-001";
    spillParams.data_file = spillBam;
    spillParams.compression = FileCompression.NONE.name();
    spillParams.files = Arrays.asList(fastqFile1.toString(), fastqFile2.toString());
    spillParams.spill_page_size = 1;
    spillParams.spill_page_size_bytes = Long.MAX_VALUE;
    spillParams.spill_abandon_limit_bytes = Long.MAX_VALUE;

    Fastq2Sam spillF2s = new Fastq2Sam();
    spillF2s.create(spillParams);

    Assert.assertEquals(8, spillF2s.getTotalReadCount());
    assertBamMatchesFastqInputs(spillParams);

    Assert.assertEquals(
        "Cross-page spill BAM should contain the same logical records as no-spill BAM",
        noSpillRecords,
        canonicalBamRecords(new File(spillBam), spillParams));
  }

  /**
   * Spill with orphans: equal-length files where some reads lack mates, combined with tiny spill
   * threshold. Verifies orphaned reads survive the spill reassembly process and the BAM is
   * identical to no-spill.
   */
  @Test
  public void testSpillReassemblyWithOrphans()
      throws IOException, ConverterException, ReadWriterException {
    Path tempDir = Files.createTempDirectory("fastq_spill_test");

    Path fastqFile1 = tempDir.resolve("input_1.fastq");
    Path fastqFile2 = tempDir.resolve("input_2.fastq");

    // Both files have 5 reads. A, C, E pair normally.
    // B and D have /1 in file1 but /2 with different keys in file2 → orphans.
    List<String> lines1 = new ArrayList<>();
    List<String> lines2 = new ArrayList<>();
    String[] names = {"READA", "READB", "READC", "READD", "READE"};
    for (String name : names) {
      lines1.addAll(
          Arrays.asList(
              "@" + name + "/1",
              "ACGTACGTACGTACGTACGTACGTACGTACGT",
              "+",
              "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"));
    }
    // File 2: A, C, E have matching mates; XORPHAN1 and XORPHAN2 are unmatched
    lines2.addAll(
        Arrays.asList(
            "@READA/2",
            "TTAATTAATTAATTAATTAATTAATTAATTAA",
            "+",
            "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"));
    lines2.addAll(
        Arrays.asList(
            "@XORPHAN1/2",
            "TTAATTAATTAATTAATTAATTAATTAATTAA",
            "+",
            "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"));
    lines2.addAll(
        Arrays.asList(
            "@READC/2",
            "TTAATTAATTAATTAATTAATTAATTAATTAA",
            "+",
            "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"));
    lines2.addAll(
        Arrays.asList(
            "@XORPHAN2/2",
            "TTAATTAATTAATTAATTAATTAATTAATTAA",
            "+",
            "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"));
    lines2.addAll(
        Arrays.asList(
            "@READE/2",
            "TTAATTAATTAATTAATTAATTAATTAATTAA",
            "+",
            "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"));

    Files.write(fastqFile1, lines1);
    Files.write(fastqFile2, lines2);

    // No-spill baseline
    // 3 pairs (A, C, E) + 4 orphans (B/1, D/1, XORPHAN1/2, XORPHAN2/2) = 10 reads
    String noSpillBam = Files.createTempFile(tempDir, "nospill", ".bam").toString();
    Fastq2Sam.Params noSpillParams = new Fastq2Sam.Params();
    noSpillParams.tmp_root = tempDir.toString();
    noSpillParams.sample_name = "SM-001";
    noSpillParams.data_file = noSpillBam;
    noSpillParams.compression = FileCompression.NONE.name();
    noSpillParams.files = Arrays.asList(fastqFile1.toString(), fastqFile2.toString());

    Fastq2Sam noSpillF2s = new Fastq2Sam();
    noSpillF2s.create(noSpillParams);

    Assert.assertEquals(10, noSpillF2s.getTotalReadCount());
    assertBamMatchesFastqInputs(noSpillParams);
    List<String> noSpillRecords = canonicalBamRecords(new File(noSpillBam), noSpillParams);

    // Spill with page_size=2
    String spillBam = Files.createTempFile(tempDir, "spill", ".bam").toString();
    Fastq2Sam.Params spillParams = new Fastq2Sam.Params();
    spillParams.tmp_root = tempDir.toString();
    spillParams.sample_name = "SM-001";
    spillParams.data_file = spillBam;
    spillParams.compression = FileCompression.NONE.name();
    spillParams.files = Arrays.asList(fastqFile1.toString(), fastqFile2.toString());
    spillParams.spill_page_size = 2;
    spillParams.spill_page_size_bytes = Long.MAX_VALUE;
    spillParams.spill_abandon_limit_bytes = Long.MAX_VALUE;

    Fastq2Sam spillF2s = new Fastq2Sam();
    spillF2s.create(spillParams);

    Assert.assertEquals(10, spillF2s.getTotalReadCount());
    assertBamMatchesFastqInputs(spillParams);

    Assert.assertEquals(
        "Spill with orphans BAM should contain the same logical records as no-spill BAM",
        noSpillRecords,
        canonicalBamRecords(new File(spillBam), spillParams));
  }

  /**
   * Large-scale spill: 20 pairs with spill_page_size=3, forcing many spill/reassembly generations.
   * Verifies all reads survive and the BAM is identical to the no-spill baseline.
   */
  @Test
  public void testSpillReassemblyManyReads()
      throws IOException, ConverterException, ReadWriterException {
    Path tempDir = Files.createTempDirectory("fastq_spill_test");

    Path fastqFile1 = tempDir.resolve("input_1.fastq");
    Path fastqFile2 = tempDir.resolve("input_2.fastq");

    int numPairs = 20;
    List<String> lines1 = new ArrayList<>();
    List<String> lines2 = new ArrayList<>();
    for (int i = 0; i < numPairs; i++) {
      String name = String.format("READ%03d", i);
      lines1.addAll(
          Arrays.asList(
              "@" + name + "/1",
              "ACGTACGTACGTACGTACGTACGTACGTACGT",
              "+",
              "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"));
      lines2.addAll(
          Arrays.asList(
              "@" + name + "/2",
              "TTAATTAATTAATTAATTAATTAATTAATTAA",
              "+",
              "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"));
    }
    Files.write(fastqFile1, lines1);
    Files.write(fastqFile2, lines2);

    // No-spill baseline
    String noSpillBam = Files.createTempFile(tempDir, "nospill", ".bam").toString();
    Fastq2Sam.Params noSpillParams = new Fastq2Sam.Params();
    noSpillParams.tmp_root = tempDir.toString();
    noSpillParams.sample_name = "SM-001";
    noSpillParams.data_file = noSpillBam;
    noSpillParams.compression = FileCompression.NONE.name();
    noSpillParams.files = Arrays.asList(fastqFile1.toString(), fastqFile2.toString());

    Fastq2Sam noSpillF2s = new Fastq2Sam();
    noSpillF2s.create(noSpillParams);

    Assert.assertEquals(numPairs * 2, noSpillF2s.getTotalReadCount());
    assertBamMatchesFastqInputs(noSpillParams);
    List<String> noSpillRecords = canonicalBamRecords(new File(noSpillBam), noSpillParams);

    // Spill with page_size=3
    String spillBam = Files.createTempFile(tempDir, "spill", ".bam").toString();
    Fastq2Sam.Params spillParams = new Fastq2Sam.Params();
    spillParams.tmp_root = tempDir.toString();
    spillParams.sample_name = "SM-001";
    spillParams.data_file = spillBam;
    spillParams.compression = FileCompression.NONE.name();
    spillParams.files = Arrays.asList(fastqFile1.toString(), fastqFile2.toString());
    spillParams.spill_page_size = 3;
    spillParams.spill_page_size_bytes = Long.MAX_VALUE;
    spillParams.spill_abandon_limit_bytes = Long.MAX_VALUE;

    Fastq2Sam spillF2s = new Fastq2Sam();
    spillF2s.create(spillParams);

    Assert.assertEquals(numPairs * 2, spillF2s.getTotalReadCount());
    assertBamMatchesFastqInputs(spillParams);

    Assert.assertEquals(
        "Many-reads spill BAM should contain the same logical records as no-spill BAM",
        noSpillRecords,
        canonicalBamRecords(new File(spillBam), spillParams));
  }

  @Test
  public void testSpillFilesUseOutputDirectoryInsteadOfTmpRoot()
      throws IOException, ConverterException, ReadWriterException {
    Path tempDir = Files.createTempDirectory("fastq_spill_location_test");
    Path spillRoot = Files.createDirectory(tempDir.resolve("spill-root"));
    Path outputDir = Files.createDirectory(tempDir.resolve("output-dir"));

    Path fastqFile1 = tempDir.resolve("input_1.fastq");
    Path fastqFile2 = tempDir.resolve("input_2.fastq");

    List<String> lines1 = new ArrayList<>();
    List<String> lines2 = new ArrayList<>();
    for (int i = 0; i < 4; i++) {
      String name = String.format("READ%03d", i);
      lines1.addAll(Arrays.asList("@" + name + "/1", "ACGTACGTACGTACGT", "+", "FFFFFFFFFFFFFFFF"));
    }
    for (int i = 3; i >= 0; i--) {
      String name = String.format("READ%03d", i);
      lines2.addAll(Arrays.asList("@" + name + "/2", "TTAATTAATTAATTAA", "+", "FFFFFFFFFFFFFFFF"));
    }
    Files.write(fastqFile1, lines1);
    Files.write(fastqFile2, lines2);

    Fastq2Sam.Params params = new Fastq2Sam.Params();
    params.tmp_root = spillRoot.toString();
    params.sample_name = "SM-001";
    params.data_file = outputDir.resolve("output.bam").toString();
    params.compression = FileCompression.NONE.name();
    params.files = Arrays.asList(fastqFile1.toString(), fastqFile2.toString());
    params.spill_page_size = 1;
    params.spill_page_size_bytes = Long.MAX_VALUE;
    params.spill_abandon_limit_bytes = Long.MAX_VALUE;

    new Fastq2Sam().create(params);

    try (Stream<Path> outputFiles = Files.list(outputDir);
        Stream<Path> spillFiles = Files.list(spillRoot)) {
      Assert.assertTrue(
          outputFiles
              .map(path -> path.getFileName().toString())
              .anyMatch(name -> name.startsWith("THREAD_") && name.contains("_PAGE_")));
      Assert.assertFalse(
          spillFiles
              .map(path -> path.getFileName().toString())
              .anyMatch(name -> name.startsWith("THREAD_") && name.contains("_PAGE_")));
    }
  }

  private Map<String, List<FastqRecord>> createFastqRecordMap(File file1, File file2) {
    Map<String, List<FastqRecord>> fastqRecordMap = new LinkedHashMap<>();

    try (FastqReader reader = new FastqReader(file1)) {
      for (FastqRecord rec : reader) {
        fastqRecordMap.put(
            rec.getReadName().substring(0, rec.getReadName().length() - 2),
            new ArrayList<>(Arrays.asList(rec)));
      }
    }

    try (FastqReader reader = new FastqReader(file2)) {
      for (FastqRecord rec : reader) {
        String key = rec.getReadName().substring(0, rec.getReadName().length() - 2);
        List<FastqRecord> val = fastqRecordMap.get(key);
        if (val == null) {
          val = new ArrayList<>(Arrays.asList(rec));
          fastqRecordMap.put(key, val);
        } else {
          // The read data in both mates should not be similar or there is no point of this test
          // otherwise.
          Assert.assertNotEquals(rec.getReadString(), val.get(0).getReadString());
          val.add(rec);
        }
      }
    }

    return fastqRecordMap;
  }

  // ---- Pipeline read name transformation tests ----
  // These tests document how the Fastq2Sam → Sam2Fastq pipeline transforms read names
  // for various pair separator conventions.

  @Test
  public void testPipelineReadNameSlashSeparator() throws Exception {
    // Input: READ/1, READ/2 — standard slash separator
    List<FastqRecord> out1 = new ArrayList<>();
    List<FastqRecord> out2 = new ArrayList<>();
    runPipelineRoundTrip(
        "@MYREAD/1\nACGTACGT\n+\nIIIIIIII\n",
        "@MYREAD/2\nTTAATTAA\n+\nKKKKKKKK\n",
        "ERR001",
        out1,
        out2);

    Assert.assertEquals(1, out1.size());
    Assert.assertEquals(1, out2.size());
    // Pipeline output: prefix.counter originalName/pairIndex
    Assert.assertEquals("ERR001.1 MYREAD/1", out1.get(0).getReadName());
    Assert.assertEquals("ERR001.1 MYREAD/2", out2.get(0).getReadName());
  }

  @Test
  public void testPipelineReadNameDotSeparator() throws Exception {
    // Input: READ.1, READ.2 — dot separator
    List<FastqRecord> out1 = new ArrayList<>();
    List<FastqRecord> out2 = new ArrayList<>();
    runPipelineRoundTrip(
        "@MYREAD.1\nACGTACGT\n+\nIIIIIIII\n",
        "@MYREAD.2\nTTAATTAA\n+\nKKKKKKKK\n",
        "ERR002",
        out1,
        out2);

    Assert.assertEquals(1, out1.size());
    Assert.assertEquals(1, out2.size());
    // Pipeline strips .N and re-emits as /N
    Assert.assertEquals("ERR002.1 MYREAD/1", out1.get(0).getReadName());
    Assert.assertEquals("ERR002.1 MYREAD/2", out2.get(0).getReadName());
  }

  @Test
  public void testPipelineReadNameColonSeparator() throws Exception {
    // Input: READ:1, READ:2 — colon separator
    List<FastqRecord> out1 = new ArrayList<>();
    List<FastqRecord> out2 = new ArrayList<>();
    runPipelineRoundTrip(
        "@MYREAD:1\nACGTACGT\n+\nIIIIIIII\n",
        "@MYREAD:2\nTTAATTAA\n+\nKKKKKKKK\n",
        "ERR003",
        out1,
        out2);

    Assert.assertEquals(1, out1.size());
    Assert.assertEquals(1, out2.size());
    Assert.assertEquals("ERR003.1 MYREAD/1", out1.get(0).getReadName());
    Assert.assertEquals("ERR003.1 MYREAD/2", out2.get(0).getReadName());
  }

  @Test
  public void testPipelineReadNameUnderscoreSeparator() throws Exception {
    // Input: READ_1, READ_2 — underscore separator
    List<FastqRecord> out1 = new ArrayList<>();
    List<FastqRecord> out2 = new ArrayList<>();
    runPipelineRoundTrip(
        "@MYREAD_1\nACGTACGT\n+\nIIIIIIII\n",
        "@MYREAD_2\nTTAATTAA\n+\nKKKKKKKK\n",
        "ERR004",
        out1,
        out2);

    Assert.assertEquals(1, out1.size());
    Assert.assertEquals(1, out2.size());
    Assert.assertEquals("ERR004.1 MYREAD/1", out1.get(0).getReadName());
    Assert.assertEquals("ERR004.1 MYREAD/2", out2.get(0).getReadName());
  }

  @Test
  public void testPipelineReadNameCasava() throws Exception {
    // Input: Casava 1.8 format — instrument:run:flowcell:lane:tile:x:y pairIndex:filter:0:barcode
    List<FastqRecord> out1 = new ArrayList<>();
    List<FastqRecord> out2 = new ArrayList<>();
    runPipelineRoundTrip(
        "@A00500:310:HG3JFDRXY:1:2101:1000:2000 1:N:0:ATCACG\nACGTACGT\n+\nIIIIIIII\n",
        "@A00500:310:HG3JFDRXY:1:2101:1000:2000 2:N:0:ATCACG\nTTAATTAA\n+\nKKKKKKKK\n",
        "ERR005",
        out1,
        out2);

    Assert.assertEquals(1, out1.size());
    Assert.assertEquals(1, out2.size());
    // Pipeline strips the Casava metadata tail (barcode/filter) through BAM round-trip,
    // keeping only the instrument part as the read name
    Assert.assertEquals(
        "ERR005.1 A00500:310:HG3JFDRXY:1:2101:1000:2000/1", out1.get(0).getReadName());
    Assert.assertEquals(
        "ERR005.1 A00500:310:HG3JFDRXY:1:2101:1000:2000/2", out2.get(0).getReadName());
  }

  @Test
  public void testPipelineReadNameNoPrefix() throws Exception {
    // No prefix — read names should just be originalName/pairIndex
    List<FastqRecord> out1 = new ArrayList<>();
    List<FastqRecord> out2 = new ArrayList<>();
    runPipelineRoundTrip(
        "@MYREAD/1\nACGTACGT\n+\nIIIIIIII\n",
        "@MYREAD/2\nTTAATTAA\n+\nKKKKKKKK\n",
        null,
        out1,
        out2);

    Assert.assertEquals(1, out1.size());
    Assert.assertEquals(1, out2.size());
    Assert.assertEquals("MYREAD/1", out1.get(0).getReadName());
    Assert.assertEquals("MYREAD/2", out2.get(0).getReadName());
  }

  @Test
  public void testPipelineReadNameMultiDigitIndex() throws Exception {
    // Input: READ.10, READ.20 — multi-digit pair index (non-standard but parseable)
    List<FastqRecord> out1 = new ArrayList<>();
    List<FastqRecord> out2 = new ArrayList<>();
    runPipelineRoundTrip(
        "@MYREAD.10\nACGTACGT\n+\nIIIIIIII\n",
        "@MYREAD.20\nTTAATTAA\n+\nKKKKKKKK\n",
        "ERR006",
        out1,
        out2);

    // Multi-digit indices: the pipeline parses 10 and 20 as pair numbers.
    // Sam2Fastq normalizes them to /1 and /2 based on the BAM first/second-of-pair flag.
    Assert.assertEquals(1, out1.size());
    Assert.assertEquals(1, out2.size());
    Assert.assertEquals("ERR006.1 MYREAD/1", out1.get(0).getReadName());
    Assert.assertEquals("ERR006.1 MYREAD/2", out2.get(0).getReadName());
  }

  /** Helper: runs Fastq2Sam → Sam2Fastq and collects output records. */
  private void runPipelineRoundTrip(
      String fastq1Content,
      String fastq2Content,
      String prefix,
      List<FastqRecord> outRecords1,
      List<FastqRecord> outRecords2)
      throws Exception {

    Path tempDir = Files.createTempDirectory("pipeline_name_test");
    Path inputFile1 = tempDir.resolve("input_1.fastq");
    Path inputFile2 = tempDir.resolve("input_2.fastq");

    Files.write(inputFile1, fastq1Content.getBytes(java.nio.charset.StandardCharsets.UTF_8));
    Files.write(inputFile2, fastq2Content.getBytes(java.nio.charset.StandardCharsets.UTF_8));

    // Fastq2Sam
    Fastq2Sam.Params f2sParams = new Fastq2Sam.Params();
    f2sParams.tmp_root = tempDir.toString();
    f2sParams.sample_name = "SM-001";
    f2sParams.data_file = tempDir.resolve("output.bam").toString();
    f2sParams.compression = FileCompression.NONE.name();
    f2sParams.files = Arrays.asList(inputFile1.toString(), inputFile2.toString());

    Fastq2Sam f2s = new Fastq2Sam();
    f2s.create(f2sParams);

    // Sam2Fastq
    Sam2Fastq.Params s2fParams = new Sam2Fastq.Params();
    s2fParams.samFile = new File(f2sParams.data_file);
    s2fParams.fastqBaseName = tempDir.resolve("result").toString();
    s2fParams.prefix = prefix;
    s2fParams.nofStreams = 3;

    Sam2Fastq s2f = new Sam2Fastq();
    s2f.create(s2fParams);

    File outFile1 = new File(s2fParams.fastqBaseName + "_1.fastq");
    File outFile2 = new File(s2fParams.fastqBaseName + "_2.fastq");

    if (outFile1.exists()) {
      try (FastqReader reader = new FastqReader(outFile1)) {
        for (FastqRecord rec : reader) outRecords1.add(rec);
      }
    }
    if (outFile2.exists()) {
      try (FastqReader reader = new FastqReader(outFile2)) {
        for (FastqRecord rec : reader) outRecords2.add(rec);
      }
    }
  }

  private static final String DEFAULT_READ_GROUP_NAME = "A";
  private static final Pattern PAIR_SUFFIX = Pattern.compile("^(.*)[.:/_]([0-9]+)$");
  private static final Pattern CASAVA_LIKE_NAME_WITHOUT_READ_INDEX =
      Pattern.compile("^([a-zA-Z0-9_-]+:[0-9]+:[a-zA-Z0-9_-]+:[0-9]+:[0-9]+:[0-9-]+:[0-9-]+)$");
  private static final Pattern CASAVA_18_NAME =
      Pattern.compile("^(.+)( +|\\t+)([0-9]+):([YN]):([0-9]*[02468])($|:.*$)");

  private static void assertBamMatchesFastqInputs(Fastq2Sam.Params params) throws IOException {
    Assert.assertEquals(
        expectedBamRecords(params), canonicalBamRecords(new File(params.data_file), params));
  }

  private static List<String> canonicalBamRecords(File bamFile, Fastq2Sam.Params params)
      throws IOException {
    List<String> records = new ArrayList<>();

    try (SamReader samReader = SamReaderFactory.makeDefault().open(bamFile)) {
      Assert.assertNotNull(samReader.getFileHeader().getReadGroup(DEFAULT_READ_GROUP_NAME));
      Assert.assertEquals(
          params.sample_name,
          samReader.getFileHeader().getReadGroup(DEFAULT_READ_GROUP_NAME).getSample());

      for (SAMRecord rec : samReader) {
        records.add(canonicalRecord(rec));
      }
    }

    Collections.sort(records);
    return records;
  }

  private static List<String> expectedBamRecords(Fastq2Sam.Params params) throws IOException {
    QualityNormalizer qualityNormalizer =
        Utils.getQualityNormalizer(
            Utils.detectFastqQualityFormat(
                params.files.get(0), params.files.size() == 2 ? params.files.get(1) : null));

    List<String> records =
        params.files.size() == 1
            ? expectedSingleEndRecords(params, qualityNormalizer)
            : expectedPairedEndRecords(params, qualityNormalizer);

    Collections.sort(records);
    return records;
  }

  private static List<String> expectedSingleEndRecords(
      Fastq2Sam.Params params, QualityNormalizer qualityNormalizer) throws IOException {
    List<String> records = new ArrayList<>();
    for (FastqInputRecord read :
        readFastqRecords(
            params.files.get(0), FileCompression.valueOf(params.compression), 1, false)) {
      records.add(
          canonicalRecord(
              read.samReadName,
              outputBases(read.bases, params.convertUracil),
              normalizedQualities(read.qualityScores, qualityNormalizer),
              false,
              false,
              false,
              true,
              false,
              read.filtered,
              DEFAULT_READ_GROUP_NAME,
              read.barcode));
    }
    return records;
  }

  private static List<String> expectedPairedEndRecords(
      Fastq2Sam.Params params, QualityNormalizer qualityNormalizer) throws IOException {
    Map<String, List<FastqInputRecord>> readsByKey = new LinkedHashMap<>();
    FileCompression compression = FileCompression.valueOf(params.compression);

    for (int fileIndex = 0; fileIndex < params.files.size(); fileIndex++) {
      for (FastqInputRecord read :
          readFastqRecords(params.files.get(fileIndex), compression, fileIndex + 1, true)) {
        readsByKey.computeIfAbsent(read.samReadName, key -> new ArrayList<>()).add(read);
      }
    }

    List<String> records = new ArrayList<>();
    for (List<FastqInputRecord> reads : readsByKey.values()) {
      Assert.assertTrue("Unexpected number of reads for spot: " + reads, reads.size() <= 2);

      if (reads.size() == 1) {
        FastqInputRecord read = reads.get(0);
        records.add(
            canonicalRecord(
                read.samReadName,
                outputBases(read.bases, params.convertUracil),
                normalizedQualities(read.qualityScores, qualityNormalizer),
                false,
                false,
                false,
                true,
                false,
                read.filtered,
                DEFAULT_READ_GROUP_NAME,
                read.barcode));
      } else {
        List<FastqInputRecord> sortedReads = new ArrayList<>(reads);
        sortedReads.sort((left, right) -> Integer.compare(left.pairNumber, right.pairNumber));

        FastqInputRecord first = sortedReads.get(0);
        records.add(
            canonicalRecord(
                first.samReadName,
                outputBases(first.bases, params.convertUracil),
                normalizedQualities(first.qualityScores, qualityNormalizer),
                true,
                true,
                false,
                true,
                true,
                first.filtered,
                DEFAULT_READ_GROUP_NAME,
                first.barcode));

        FastqInputRecord second = sortedReads.get(1);
        records.add(
            canonicalRecord(
                second.samReadName,
                outputBases(second.bases, params.convertUracil),
                normalizedQualities(second.qualityScores, qualityNormalizer),
                true,
                false,
                true,
                true,
                true,
                second.filtered,
                DEFAULT_READ_GROUP_NAME,
                second.barcode));
      }
    }
    return records;
  }

  private static List<FastqInputRecord> readFastqRecords(
      String fastqFile, FileCompression compression, int defaultReadIndex, boolean paired)
      throws IOException {
    List<FastqInputRecord> records = new ArrayList<>();

    try (InputStream is = compression.open(fastqFile, false);
        BufferedReader reader =
            new BufferedReader(new InputStreamReader(is, StandardCharsets.UTF_8))) {
      String nameLine;
      while ((nameLine = readNextNonBlankLine(reader)) != null) {
        String bases = reader.readLine();
        String separator = reader.readLine();
        String qualityScores = reader.readLine();

        Assert.assertTrue(
            "FASTQ read name line should start with @: " + nameLine, nameLine.startsWith("@"));
        Assert.assertNotNull("Missing FASTQ bases line after: " + nameLine, bases);
        Assert.assertNotNull("Missing FASTQ separator line after: " + nameLine, separator);
        Assert.assertNotNull("Missing FASTQ quality line after: " + nameLine, qualityScores);
        Assert.assertTrue(
            "FASTQ separator line should start with +: " + separator, separator.startsWith("+"));

        String readName = nameLine.substring(1);
        ReadNameInfo readNameInfo = readNameInfo(readName, defaultReadIndex, paired);
        records.add(
            new FastqInputRecord(
                readName,
                readNameInfo.samReadName,
                readNameInfo.pairNumber,
                bases,
                qualityScores,
                readNameInfo.filtered,
                readNameInfo.barcode));
      }
    }

    return records;
  }

  private static String readNextNonBlankLine(BufferedReader reader) throws IOException {
    String line;
    while ((line = reader.readLine()) != null) {
      if (!line.trim().isEmpty()) {
        return line;
      }
    }
    return null;
  }

  private static ReadNameInfo readNameInfo(String readName, int defaultReadIndex, boolean paired) {
    Matcher casava18Matcher = CASAVA_18_NAME.matcher(readName);
    if (casava18Matcher.matches()) {
      String barcode =
          casava18Matcher.group(6).isEmpty() ? null : casava18Matcher.group(6).substring(1);
      return new ReadNameInfo(
          casava18Matcher.group(1),
          Integer.parseInt(casava18Matcher.group(3)),
          "Y".equals(casava18Matcher.group(4)),
          barcode == null || barcode.isEmpty() ? null : barcode);
    }

    if (!paired) {
      int slashIndex = readName.lastIndexOf('/');
      return new ReadNameInfo(
          slashIndex == -1 ? readName : readName.substring(0, slashIndex),
          defaultReadIndex,
          false,
          null);
    }

    if (!CASAVA_LIKE_NAME_WITHOUT_READ_INDEX.matcher(readName).matches()) {
      Matcher pairSuffixMatcher = PAIR_SUFFIX.matcher(readName);
      if (pairSuffixMatcher.matches()) {
        return new ReadNameInfo(
            pairSuffixMatcher.group(1), Integer.parseInt(pairSuffixMatcher.group(2)), false, null);
      }
    }

    return new ReadNameInfo(readName, defaultReadIndex, false, null);
  }

  private static String outputBases(String bases, boolean convertUracil) {
    return (convertUracil ? Utils.replaceUracilBases(bases) : bases)
        .toUpperCase(java.util.Locale.ROOT);
  }

  private static byte[] normalizedQualities(
      String qualityScores, QualityNormalizer qualityNormalizer) {
    byte[] normalizedQualities = qualityScores.getBytes(StandardCharsets.UTF_8);
    qualityNormalizer.normalize(normalizedQualities);
    return normalizedQualities;
  }

  private static String canonicalRecord(SAMRecord record) {
    boolean paired = record.getReadPairedFlag();
    return canonicalRecord(
        record.getReadName(),
        record.getReadString(),
        record.getBaseQualities(),
        paired,
        paired && record.getFirstOfPairFlag(),
        paired && record.getSecondOfPairFlag(),
        record.getReadUnmappedFlag(),
        paired && record.getMateUnmappedFlag(),
        record.getReadFailsVendorQualityCheckFlag(),
        String.valueOf(record.getAttribute("RG")),
        String.valueOf(record.getAttribute("BC")));
  }

  private static String canonicalRecord(
      String readName,
      String bases,
      byte[] baseQualities,
      boolean paired,
      boolean firstOfPair,
      boolean secondOfPair,
      boolean unmapped,
      boolean mateUnmapped,
      boolean qcFail,
      String readGroup,
      String barcode) {
    return String.join(
        "\t",
        "name=" + readName,
        "bases=" + bases,
        "qualities=" + Arrays.toString(baseQualities),
        "paired=" + paired,
        "first=" + firstOfPair,
        "second=" + secondOfPair,
        "unmapped=" + unmapped,
        "mateUnmapped=" + mateUnmapped,
        "qcFail=" + qcFail,
        "RG=" + readGroup,
        "BC=" + barcode);
  }

  private static final class ReadNameInfo {
    final String samReadName;
    final int pairNumber;
    final boolean filtered;
    final String barcode;

    ReadNameInfo(String samReadName, int pairNumber, boolean filtered, String barcode) {
      this.samReadName = samReadName;
      this.pairNumber = pairNumber;
      this.filtered = filtered;
      this.barcode = barcode;
    }
  }

  private static final class FastqInputRecord {
    final String originalReadName;
    final String samReadName;
    final int pairNumber;
    final String bases;
    final String qualityScores;
    final boolean filtered;
    final String barcode;

    FastqInputRecord(
        String originalReadName,
        String samReadName,
        int pairNumber,
        String bases,
        String qualityScores,
        boolean filtered,
        String barcode) {
      this.originalReadName = originalReadName;
      this.samReadName = samReadName;
      this.pairNumber = pairNumber;
      this.bases = bases;
      this.qualityScores = qualityScores;
      this.filtered = filtered;
      this.barcode = barcode;
    }

    @Override
    public String toString() {
      return originalReadName;
    }
  }
}
