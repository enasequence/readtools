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
package uk.ac.ebi.ena.readtools.v2;

import static org.junit.Assert.*;
import static uk.ac.ebi.ena.readtools.v2.TestFileUtil.READ_COUNT_LIMIT;
import static uk.ac.ebi.ena.readtools.v2.TestFileUtil.createOutputFolder;
import static uk.ac.ebi.ena.readtools.v2.validator.InsdcReadsValidator.ERROR_BASES_QUALITIES_LENGTH_MISMATCH;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import org.junit.Test;
import uk.ac.ebi.ena.readtools.v2.provider.ReadsProviderFactory;
import uk.ac.ebi.ena.readtools.v2.validator.ReadsValidationException;
import uk.ac.ebi.ena.readtools.v2.validator.ValidatorWrapper;

public class ValidatorWrapperTest {

  // Verifies BAM wrapper flow validates and collects stats for a single valid file.
  @Test
  public void collectsStatsForSingleValidBamFile() throws IOException, ReadsValidationException {
    File outputDir = createOutputFolder();
    Path bamPath = outputDir.toPath().resolve("single-valid.bam");
    writeBamWithRecord(bamPath, createUnmappedRecord("read-1", "ACGT", "!!!!"));

    ValidatorWrapper wrapper =
        new ValidatorWrapper(
            Collections.singletonList(bamPath.toFile()), FileFormat.BAM, READ_COUNT_LIMIT);
    wrapper.run();

    List<ValidatorWrapper.FileQualityStats> fileQualityStats = wrapper.getFileQualityStats();
    assertEquals(1, fileQualityStats.size());
    assertEquals(bamPath.toFile(), fileQualityStats.get(0).getFile());
    assertEquals(1, fileQualityStats.get(0).getReadCount());
  }

  // Verifies BAM wrapper flow validates and collects stats for multiple files.
  @Test
  public void collectsStatsForMultipleValidBamFiles() throws IOException, ReadsValidationException {
    File outputDir = createOutputFolder();
    Path bamPath1 = outputDir.toPath().resolve("multi-valid-1.bam");
    Path bamPath2 = outputDir.toPath().resolve("multi-valid-2.bam");
    writeBamWithRecord(bamPath1, createUnmappedRecord("read-1", "ACGT", "!!!!"));
    writeBamWithRecord(bamPath2, createUnmappedRecord("read-2", "TGCA", "!!!!"));

    ValidatorWrapper wrapper =
        new ValidatorWrapper(
            Arrays.asList(bamPath1.toFile(), bamPath2.toFile()), FileFormat.BAM, READ_COUNT_LIMIT);
    wrapper.run();

    List<ValidatorWrapper.FileQualityStats> fileQualityStats = wrapper.getFileQualityStats();
    assertEquals(2, fileQualityStats.size());
    assertEquals(bamPath1.toFile(), fileQualityStats.get(0).getFile());
    assertEquals(bamPath2.toFile(), fileQualityStats.get(1).getFile());
    assertEquals(1, fileQualityStats.get(0).getReadCount());
    assertEquals(1, fileQualityStats.get(1).getReadCount());
  }

  // Verifies CRAM wrapper flow validates and collects stats for a valid file.
  @Test
  public void collectsStatsForSingleValidCramFile() throws IOException, ReadsValidationException {
    File outputDir = createOutputFolder();
    Path cramPath = outputDir.toPath().resolve("single-valid.cram");
    writeCramWithRecord(cramPath, createUnmappedRecord("read-1", "ACGT", "!!!!"));

    ValidatorWrapper wrapper =
        new ValidatorWrapper(
            Collections.singletonList(cramPath.toFile()), FileFormat.CRAM, READ_COUNT_LIMIT);
    wrapper.run();

    List<ValidatorWrapper.FileQualityStats> fileQualityStats = wrapper.getFileQualityStats();
    assertEquals(1, fileQualityStats.size());
    assertEquals(cramPath.toFile(), fileQualityStats.get(0).getFile());
    assertEquals(1, fileQualityStats.get(0).getReadCount());
  }

  // Verifies BAM wrapper flow propagates validation errors and keeps stats from files validated
  // earlier.
  @Test
  public void propagatesAndPreservesPartialStatsWhenBamValidationFailsMidRun() throws IOException {
    File outputDir = createOutputFolder();
    Path validBamPath = outputDir.toPath().resolve("valid-first.bam");
    Path invalidBamPath = outputDir.toPath().resolve("invalid-second.bam");
    writeBamWithRecord(validBamPath, createUnmappedRecord("read-valid", "ACGT", "!!!!"));
    writeBamWithRecord(
        invalidBamPath, createUnmappedRecordWithoutQualities("read-invalid", "ACGT"));

    ValidatorWrapper wrapper =
        new ValidatorWrapper(
            Arrays.asList(validBamPath.toFile(), invalidBamPath.toFile()),
            FileFormat.BAM,
            READ_COUNT_LIMIT);

    try {
      wrapper.run();
      fail();
    } catch (ReadsValidationException validationException) {
      assertTrue(
          validationException.getErrorMessage().contains(ERROR_BASES_QUALITIES_LENGTH_MISMATCH));
    }

    List<ValidatorWrapper.FileQualityStats> fileQualityStats = wrapper.getFileQualityStats();
    assertEquals(1, fileQualityStats.size());
    assertEquals(validBamPath.toFile(), fileQualityStats.get(0).getFile());
    assertEquals(1, fileQualityStats.get(0).getReadCount());
  }

  // Verifies FASTA is explicitly unsupported in ReadsProviderFactory.
  @Test
  public void readsProviderFactoryRejectsFastaFormatAsNotImplemented() {
    ReadsProviderFactory readsProviderFactory =
        new ReadsProviderFactory(new File("unused"), FileFormat.FASTA);

    try {
      readsProviderFactory.makeReadsProvider();
      fail();
    } catch (ReadsValidationException validationException) {
      assertTrue(validationException.getMessage().contains("not implemented"));
    }
  }

  // Verifies FASTA is explicitly unsupported in ValidatorWrapper.
  @Test
  public void validatorWrapperRejectsFastaFormatAsNotImplemented() {
    ValidatorWrapper wrapper =
        new ValidatorWrapper(Collections.emptyList(), FileFormat.FASTA, READ_COUNT_LIMIT);

    try {
      wrapper.run();
      fail();
    } catch (ReadsValidationException validationException) {
      assertTrue(validationException.getMessage().contains("not implemented"));
    }
  }

  private void writeBamWithRecord(Path bamPath, SAMRecord record) {
    try (SAMFileWriter writer =
        new SAMFileWriterFactory().makeBAMWriter(record.getHeader(), false, bamPath.toFile())) {
      writer.addAlignment(record);
    }
  }

  private void writeCramWithRecord(Path cramPath, SAMRecord record) throws IOException {
    try (OutputStream outputStream = Files.newOutputStream(cramPath);
        SAMFileWriter writer =
            new SAMFileWriterFactory()
                .makeCRAMWriter(record.getHeader(), outputStream, (Path) null)) {
      writer.addAlignment(record);
    }
  }

  private SAMRecord createUnmappedRecord(String readName, String readString, String qualityString) {
    SAMFileHeader header = new SAMFileHeader();
    header.setSortOrder(SAMFileHeader.SortOrder.unsorted);

    SAMRecord record = new SAMRecord(header);
    record.setReadName(readName);
    record.setReadString(readString);
    record.setBaseQualityString(qualityString);
    record.setReadUnmappedFlag(true);
    return record;
  }

  private SAMRecord createUnmappedRecordWithoutQualities(String readName, String readString) {
    SAMRecord record = createUnmappedRecord(readName, readString, "!!!!");
    record.setBaseQualities(SAMRecord.NULL_QUALS);
    return record;
  }
}
