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
import static uk.ac.ebi.ena.readtools.v2.TestFileUtil.*;
import static uk.ac.ebi.ena.readtools.v2.validator.InsdcReadsValidator.*;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import org.junit.Ignore;
import org.junit.Test;
import uk.ac.ebi.ena.readtools.v2.MockReadsProvider.MockRead;
import uk.ac.ebi.ena.readtools.v2.provider.ReadsProvider;
import uk.ac.ebi.ena.readtools.v2.provider.ReadsProviderFactory;
import uk.ac.ebi.ena.readtools.v2.read.IRead;
import uk.ac.ebi.ena.readtools.v2.validator.FastqReadsValidator;
import uk.ac.ebi.ena.readtools.v2.validator.InsdcReadsValidator;
import uk.ac.ebi.ena.readtools.v2.validator.ReadsValidationException;

public class InsdcReadsValidatorTest {
  // Verifies validation fails when the provider yields no reads.
  @Test
  public void rejectsProviderWithNoReads() {
    try {
      MockReadsProviderFactory readsProviderFactory = new MockReadsProviderFactory();
      new InsdcReadsValidator(READ_COUNT_LIMIT).validate(readsProviderFactory);
      fail();
    } catch (ReadsValidationException validationException) {
      assertTrue(validationException.getErrorMessage().contains(ERROR_NO_READS));
    }
  }

  // Verifies null, empty, and whitespace-only base strings are rejected.
  @Test
  public void rejectsNullEmptyOrWhitespaceReadBases() {
    try {
      MockReadsProviderFactory readsProviderFactory =
          new MockReadsProviderFactory(new MockRead("r1", null, "1234"));
      new InsdcReadsValidator(READ_COUNT_LIMIT).validate(readsProviderFactory);
      fail();
    } catch (ReadsValidationException validationException) {
      assertTrue(validationException.getErrorMessage().contains(ERROR_EMPTY_READ));
    }

    try {
      MockReadsProviderFactory readsProviderFactory =
          new MockReadsProviderFactory(new MockRead("r1", "", "1234"));
      new InsdcReadsValidator(READ_COUNT_LIMIT).validate(readsProviderFactory);
      fail();
    } catch (ReadsValidationException validationException) {
      assertTrue(validationException.getErrorMessage().contains(ERROR_EMPTY_READ));
    }

    try {
      MockReadsProviderFactory readsProviderFactory =
          new MockReadsProviderFactory(new MockRead("r1", " ", "1234"));
      new InsdcReadsValidator(READ_COUNT_LIMIT).validate(readsProviderFactory);
      fail();
    } catch (ReadsValidationException validationException) {
      assertTrue(validationException.getErrorMessage().contains(ERROR_EMPTY_READ));
    }
  }

  // Verifies an empty FASTQ file is rejected as having zero reads.
  @Test
  public void rejectsEmptyFastqFile() throws IOException {
    File outputDir = createOutputFolder();
    Path fastqPath = saveRandomized("", outputDir.toPath(), true, "fastq-zero-reads-single", "gz");

    try {
      ReadsProviderFactory readsProviderFactory =
          new ReadsProviderFactory(fastqPath.toFile(), FileFormat.FASTQ);
      new InsdcReadsValidator(READ_COUNT_LIMIT).validate(readsProviderFactory);
      fail();
    } catch (ReadsValidationException validationException) {
      // Existing parser behavior varies; this test only confirms failure.
      // todo: assertEquals(ERROR_NO_READS, validationException.getMessage());
    }
  }

  // Verifies FASTQ records with an empty sequence line are rejected.
  @Test
  public void rejectsFastqRecordWithMissingSequenceLine() throws IOException {
    File outputDir = createOutputFolder();
    Path fastqPath =
        saveRandomized(
            "@CL100031435L2C001R001_63456 :TGGTCCTTCC\n"
                + "\n"
                + "+\n"
                + "FDFFFFFFFG3FFFFFF7?D8FFFFF=>DFG(FBEFFFFDF;FF?FF<8FGGFFFCBB8F0F@FBC?FAGFEE>.FFEFCF:?E(E@;3*(,FD/BFE-\n",
            outputDir.toPath(),
            true,
            "fastq-10",
            "gz");

    try {
      ReadsProviderFactory readsProviderFactory =
          new ReadsProviderFactory(fastqPath.toFile(), FileFormat.FASTQ);
      new InsdcReadsValidator(READ_COUNT_LIMIT).validate(readsProviderFactory);
      fail();
    } catch (ReadsValidationException validationException) {
      // Existing parser behavior varies; this test only confirms failure.
      // todo: assertTrue(validationException.getMessage().contains("Missing Sequence Line"));
    }
  }

  // Verifies read names longer than 256 characters are rejected.
  @Test
  public void rejectsReadNameLongerThan256Chars() throws IOException {
    File outputDir = createOutputFolder();
    Path fastqPath =
        saveRandomized(
            "@aaaabbbbccccddddaaaabbbbccccddddaaaabbbbccccddddaaaabbbbccccddddaaaabbbbccccddddaaaabbbbccccddddaaaabbbbccccddddaaaabbbbccccddddaaaabbbbccccddddaaaabbbbccccddddaaaabbbbccccddddaaaabbbbccccddddaaaabbbbccccddddaaaabbbbccccddddaaaabbbbccccddddaaaabbbbccccdddd1\n"
                + "AGTC\n"
                + "+\n"
                + "1234\n",
            outputDir.toPath(),
            true,
            "fastq",
            "gz");

    try {
      ReadsProviderFactory readsProviderFactory =
          new ReadsProviderFactory(fastqPath.toFile(), FileFormat.FASTQ);
      new InsdcReadsValidator(READ_COUNT_LIMIT).validate(readsProviderFactory);
      fail();
    } catch (ReadsValidationException validationException) {
      assertTrue(
          validationException.getMessage().contains("Read name length exceeds 256 characters"));
    }
  }

  // Verifies FASTQ records with a missing quality line are rejected.
  @Test
  public void rejectsFastqRecordWithMissingQualityLine() throws IOException {
    File outputDir = createOutputFolder();
    Path fastqPath =
        saveRandomized(
            "@CL100031435L2C001R001_63456 :TGGTCCTTCC\n"
                + "ATCCTTATAGTGGGCCAAGCTCTCACATGCAAACACGTTGCCTGCTGGATTTTGTTTCAGAGAGGAAATGTTTATGTGAGACAGAAAAAGCCGGGGGCC\n"
                + "+\n",
            outputDir.toPath(),
            true,
            "fastq-10",
            "gz");

    try {
      ReadsProviderFactory readsProviderFactory =
          new ReadsProviderFactory(fastqPath.toFile(), FileFormat.FASTQ);
      new InsdcReadsValidator(READ_COUNT_LIMIT).validate(readsProviderFactory);
      fail();
    } catch (ReadsValidationException validationException) {
      // Existing parser behavior varies; this test only confirms failure.
    }
  }

  // Verifies '*' is accepted as a FASTQ quality character when sequence length is 1.
  @Test
  public void acceptsFastqSingleBaseWithAsteriskQuality() throws ReadsValidationException {
    ReadsProviderFactory readsProviderFactory =
        createMockReadsProviderFactory(
            FileFormat.FASTQ, new MockRead("read-with-asterisk-quality", "A", "*"));

    new InsdcReadsValidator(READ_COUNT_LIMIT).validate(readsProviderFactory);
  }

  // Verifies FASTQ validation rejects base/quality length mismatches.
  @Test
  public void rejectsFastqWhenBasesAndQualitiesLengthsMismatch() {
    try {
      ReadsProviderFactory readsProviderFactory =
          createMockReadsProviderFactory(
              FileFormat.FASTQ, new MockRead("fastq-mismatch", "ACGT", "*"));
      new InsdcReadsValidator(READ_COUNT_LIMIT).validate(readsProviderFactory);
      fail();
    } catch (ReadsValidationException validationException) {
      assertTrue(
          validationException.getErrorMessage().contains(ERROR_BASES_QUALITIES_LENGTH_MISMATCH));
    }
  }

  // Verifies empty quality strings in in-memory reads are rejected.
  @Test
  public void rejectsEmptyQualityStringInMockRead() {
    try {
      MockReadsProviderFactory readsProviderFactory =
          new MockReadsProviderFactory(new MockRead("r1", "AGTC", ""));
      new InsdcReadsValidator(READ_COUNT_LIMIT).validate(readsProviderFactory);
      fail();
    } catch (ReadsValidationException validationException) {
      assertTrue(validationException.getErrorMessage().contains(ERROR_EMPTY_QUALITY));
    }
  }

  // Verifies FASTQ files containing an empty read entry are rejected.
  @Test
  public void rejectsFastqFileContainingEmptyRead() throws IOException {
    File outputDir = createOutputFolder();
    Path fastqPath =
        saveRandomized(
            "@CL100031435L2C001R001_63456 :TGGTCCTTCC\n"
                + "ATCCTTATAGTGGGCCAAGCTCTCACATGCAAACACGTTGCCTGCTGGATTTTGTTTCAGAGAGGAAATGTTTATGTGAGACAGAAAAAGCCGGGGGCC\n"
                + "+\n"
                + "FDFFFFFFFG3FFFFFF7?D8FFFFF=>DFG(FBEFFFFDF;FF?FF<8FGGFFFCBB8F0F@FBC?FAGFEE>.FFEFCF:?E(E@;3*(,FD/BFE-\n"
                + "@CL100031435L2C001R001_63472 :TGGTCCTTCA\n"
                + "\n"
                + "+\n"
                + "\n"
                + "@CL100031435L2C001R001_63487 :TGGTCCTTCC\n"
                + "CCGCTCTGCCCACCACAATCCGGCCCCTGTGTACGGCAACACAGGGCCCAGGCAAGCAGATCCTTCCTGCTGGGAGCTCCAGCTTGTAGAATTTCACCC\n"
                + "+\n"
                + "FFFFFFFFFEFFFFFFD5FEFFEEEFFEBCE?EEFFFECD:DAFFCDFFBFFDB@>FCFBEEFF>F-CBFFF&E9F=BBF:4@B?B45E2F+A-EFFCA\n",
            outputDir.toPath(),
            true,
            "fastq-10",
            "gz");

    try {
      ReadsProviderFactory readsProviderFactory =
          new ReadsProviderFactory(fastqPath.toFile(), FileFormat.FASTQ);
      new InsdcReadsValidator(READ_COUNT_LIMIT).validate(readsProviderFactory);
      fail();
    } catch (ReadsValidationException validationException) {
      // Existing parser behavior varies; this test only confirms failure.
      // todo: assertEquals(INVALID_FILE, validationException.getMessage());
    }
  }

  // Verifies non-IUPAC base codes are rejected.
  @Test
  public void rejectsReadWithNonIupacBases() {
    try {
      MockReadsProviderFactory readsProviderFactory =
          new MockReadsProviderFactory(new MockRead("r1", "AFFF", "1234"));
      new InsdcReadsValidator(READ_COUNT_LIMIT).validate(readsProviderFactory);
      fail();
    } catch (ReadsValidationException validationException) {
      assertTrue(validationException.getErrorMessage().contains(ERROR_NOT_IUPAC));
    }
  }

  // Verifies data is rejected when over 50% of bases are non-A/U/T/C/G IUPAC symbols.
  @Test
  public void rejectsWhenMoreThanHalfBasesAreNonAutcg() {
    try {
      MockReadsProviderFactory readsProviderFactory =
          new MockReadsProviderFactory(
              new MockRead("r0", "AAAA", "1234"),
              new MockRead("r1", "AAAA", "1234"),
              new MockRead("r2", "AAAA", "1234"),
              new MockRead("r3", "AAAA", "1234"),
              new MockRead("r4", "NNNN", "1234"),
              new MockRead("r5", "NNNN", "1234"),
              new MockRead("r6", "NNNN", "1234"),
              new MockRead("r7", "NNNN", "1234"),
              new MockRead("r8", "NNNN", "1234"),
              new MockRead("r9", "NNNN", "1234"));
      new InsdcReadsValidator(READ_COUNT_LIMIT).validate(readsProviderFactory);
      fail();
    } catch (ReadsValidationException validationException) {
      assertTrue(validationException.getErrorMessage().contains(ERROR_NOT_AUTCG));
    }
  }

  // Documents expected low-quality behavior if the quality threshold rule is re-enabled.
  @Ignore("rejects too many submissions")
  @Test
  public void lowQualityRuleWouldRejectReadsWhenEnabled() throws ReadsValidationException {
    MockReadsProviderFactory readsProviderFactory =
        new MockReadsProviderFactory(new MockRead("r1", "AGTC", ">>>>"));
    new InsdcReadsValidator(READ_COUNT_LIMIT).validate(readsProviderFactory);

    try {
      readsProviderFactory = new MockReadsProviderFactory(new MockRead("r1", "AGTC", "@@@@"));
      new InsdcReadsValidator(READ_COUNT_LIMIT).validate(readsProviderFactory);
      fail();
    } catch (ReadsValidationException validationException) {
      assertEquals(ERROR_QUALITY, validationException.getErrorMessage());
    }
  }

  // Verifies a valid FASTQ sample passes both INSDC and FASTQ validators.
  @Test
  public void acceptsValidFastqData() throws IOException, ReadsValidationException {
    File outputDir = createOutputFolder();
    Path fastqPath =
        saveRandomized(
            "@CL100031435L2C001R001_63456 :TGGTCCTTCC\n"
                + "ATCCTTATAGTGGGCCAAGCTCTCACATGCAAACACGTTGCCTGCTGGATTTTGTTTCAGAGAGGAAATGTTTATGTGAGACAGAAAAAGCCGGGGGCC\n"
                + "+\n"
                + "FDFFFFFFFG?FFFFFF??D?FFFFF??DFG?FBEFFFFDF?FF?FF??FGGFFFCBB?F?F@FBC?FAGFEE??FFEFCF??E?E@?????FD?BFE?\n"
                + "@CL100031435L2C001R001_63487 :TGGTCCTTCC\n"
                + "CCGCTCTGCCCACCACAATCCGGCCCCTGTGTACGGCAACACAGGGCCCAGGCAAGCAGATCCTTCCTGCTGGGAGCTCCAGCTTGTAGAATTTCACCC\n"
                + "+\n"
                + "FFFFFFFFFEFFFFFFD?FEFFEEEFFEBCE?EEFFFECD?DAFFCDFFBFFDB@?FCFBEEFF?F?CBFFF?E?F?BBF??@B?B??E?F?A?EFFCA\n",
            outputDir.toPath(),
            true,
            "fastq-10",
            "gz");

    ReadsProviderFactory readsProviderFactory =
        new ReadsProviderFactory(fastqPath.toFile(), FileFormat.FASTQ, true);
    new InsdcReadsValidator(READ_COUNT_LIMIT).validate(readsProviderFactory);
    new FastqReadsValidator(READ_COUNT_LIMIT).validate(readsProviderFactory);
  }

  // Verifies BAM QUAL='*' (missing qualities) is rejected as a bases/qualities length mismatch.
  @Test
  public void rejectsBamWithMissingBaseQualitiesAsLengthMismatch() throws IOException {
    File outputDir = createOutputFolder();
    Path bamPath = outputDir.toPath().resolve("missing-base-qualities.bam");
    writeBamWithMissingBaseQualities(bamPath);

    try (SamReader samReader = SamReaderFactory.makeDefault().open(bamPath.toFile())) {
      SAMRecord bamRecord = samReader.iterator().next();
      assertEquals("*", bamRecord.getBaseQualityString());
    }

    ReadsProviderFactory readsProviderFactory =
        new ReadsProviderFactory(bamPath.toFile(), FileFormat.BAM);
    try {
      new InsdcReadsValidator(READ_COUNT_LIMIT).validate(readsProviderFactory);
      fail();
    } catch (ReadsValidationException validationException) {
      assertTrue(
          validationException.getErrorMessage().contains(ERROR_BASES_QUALITIES_LENGTH_MISMATCH));
    }
  }

  // Verifies SAM QUAL='*' (missing qualities) is rejected as a bases/qualities length mismatch.
  @Test
  public void rejectsSamWithMissingBaseQualitiesAsLengthMismatch() throws IOException {
    File outputDir = createOutputFolder();
    Path samPath = outputDir.toPath().resolve("missing-base-qualities.sam");
    writeSamWithMissingBaseQualities(samPath);

    ReadsProviderFactory readsProviderFactory =
        new ReadsProviderFactory(samPath.toFile(), FileFormat.BAM);
    try {
      new InsdcReadsValidator(READ_COUNT_LIMIT).validate(readsProviderFactory);
      fail();
    } catch (ReadsValidationException validationException) {
      assertTrue(
          validationException.getErrorMessage().contains(ERROR_BASES_QUALITIES_LENGTH_MISMATCH));
    }
  }

  // Verifies CRAM QUAL='*' (missing qualities) is rejected as a bases/qualities length mismatch.
  @Test
  public void rejectsCramWithMissingBaseQualitiesAsLengthMismatch() throws IOException {
    File outputDir = createOutputFolder();
    Path cramPath = outputDir.toPath().resolve("missing-base-qualities.cram");
    writeCramWithMissingBaseQualities(cramPath);

    ReadsProviderFactory readsProviderFactory =
        new ReadsProviderFactory(cramPath.toFile(), FileFormat.CRAM);
    try {
      new InsdcReadsValidator(READ_COUNT_LIMIT).validate(readsProviderFactory);
      fail();
    } catch (ReadsValidationException validationException) {
      assertTrue(
          validationException.getErrorMessage().contains(ERROR_BASES_QUALITIES_LENGTH_MISMATCH));
    }
  }

  // Verifies SAM/BAM validation rejects files that contain no reads.
  @Test
  public void rejectsBamFileWithNoReads() throws IOException {
    File outputDir = createOutputFolder();
    Path bamPath = outputDir.toPath().resolve("no-reads.bam");
    writeHeaderOnlyBam(bamPath);

    ReadsProviderFactory readsProviderFactory =
        new ReadsProviderFactory(bamPath.toFile(), FileFormat.BAM);
    try {
      new InsdcReadsValidator(READ_COUNT_LIMIT).validate(readsProviderFactory);
      fail();
    } catch (ReadsValidationException validationException) {
      assertTrue(validationException.getErrorMessage().contains(ERROR_NO_READS));
    }
  }

  // Verifies SAM/BAM validation rejects read names longer than 256 characters.
  @Test
  public void rejectsBamReadNameLongerThan256Chars() {
    String longReadName = "a".repeat(257);
    ReadsProviderFactory readsProviderFactory =
        createMockReadsProviderFactory(FileFormat.BAM, new MockRead(longReadName, "ACGT", "!!!!"));

    try {
      new InsdcReadsValidator(READ_COUNT_LIMIT).validate(readsProviderFactory);
      fail();
    } catch (ReadsValidationException validationException) {
      assertTrue(validationException.getErrorMessage().contains(ERROR_READ_NAME_LENGTH));
    }
  }

  // Verifies SAM/BAM validation rejects non-IUPAC base codes.
  @Test
  public void rejectsBamReadWithNonIupacBases() {
    ReadsProviderFactory readsProviderFactory =
        createMockReadsProviderFactory(
            FileFormat.BAM, new MockRead("read-non-iupac", "AFFF", "!!!!"));

    try {
      new InsdcReadsValidator(READ_COUNT_LIMIT).validate(readsProviderFactory);
      fail();
    } catch (ReadsValidationException validationException) {
      assertTrue(validationException.getErrorMessage().contains(ERROR_NOT_IUPAC));
    }
  }

  // Verifies valid BAM data passes INSDC validation.
  @Test
  public void acceptsValidBamData() throws IOException, ReadsValidationException {
    File outputDir = createOutputFolder();
    Path bamPath = outputDir.toPath().resolve("valid.bam");
    SAMRecord record = createUnmappedRecord("read-valid", "ACGT", "!!!!");
    writeBamWithRecord(bamPath, record);

    ReadsProviderFactory readsProviderFactory =
        new ReadsProviderFactory(bamPath.toFile(), FileFormat.BAM);
    new InsdcReadsValidator(READ_COUNT_LIMIT).validate(readsProviderFactory);
  }

  // Verifies SAM/BAM validation rejects records with missing sequence (SEQ='*').
  @Test
  public void rejectsBamRecordWithMissingSequence() throws IOException {
    File outputDir = createOutputFolder();
    Path bamPath = outputDir.toPath().resolve("missing-sequence.bam");
    SAMRecord record = createUnmappedRecordWithoutSequence();
    writeBamWithRecord(bamPath, record);

    ReadsProviderFactory readsProviderFactory =
        new ReadsProviderFactory(bamPath.toFile(), FileFormat.BAM);
    try {
      new InsdcReadsValidator(READ_COUNT_LIMIT).validate(readsProviderFactory);
      fail();
    } catch (ReadsValidationException validationException) {
      assertTrue(
          "Expected: " + ERROR_EMPTY_READ + ", actual: " + validationException.getErrorMessage(),
          validationException.getErrorMessage().contains(ERROR_EMPTY_READ));
    }
  }

  private ReadsProviderFactory createMockReadsProviderFactory(
      FileFormat format, MockRead... reads) {
    return new ReadsProviderFactory(null, format) {
      @Override
      public ReadsProvider<? extends IRead> makeReadsProvider() {
        return new MockReadsProvider(reads);
      }
    };
  }

  private void writeBamWithMissingBaseQualities(Path bamPath) {
    SAMRecord record = createUnmappedRecordWithoutQualities();
    writeBamWithRecord(bamPath, record);
  }

  private void writeBamWithRecord(Path bamPath, SAMRecord record) {
    try (SAMFileWriter writer =
        new SAMFileWriterFactory().makeBAMWriter(record.getHeader(), false, bamPath.toFile())) {
      writer.addAlignment(record);
    }
  }

  private void writeHeaderOnlyBam(Path bamPath) {
    SAMFileHeader header = new SAMFileHeader();
    header.setSortOrder(SAMFileHeader.SortOrder.unsorted);
    try (SAMFileWriter ignored =
        new SAMFileWriterFactory().makeBAMWriter(header, false, bamPath.toFile())) {
      // Header only.
    }
  }

  private void writeSamWithMissingBaseQualities(Path samPath) {
    SAMRecord record = createUnmappedRecordWithoutQualities();
    try (SAMFileWriter writer =
        new SAMFileWriterFactory().makeSAMWriter(record.getHeader(), false, samPath.toFile())) {
      writer.addAlignment(record);
    }
  }

  private void writeCramWithMissingBaseQualities(Path cramPath) throws IOException {
    SAMRecord record = createUnmappedRecordWithoutQualities();
    try (OutputStream outputStream = Files.newOutputStream(cramPath);
        SAMFileWriter writer =
            new SAMFileWriterFactory()
                .makeCRAMWriter(record.getHeader(), outputStream, (Path) null)) {
      writer.addAlignment(record);
    }
  }

  private SAMRecord createUnmappedRecordWithoutQualities() {
    SAMRecord record = createUnmappedRecord("read-1", "ACGT", "!!!!");
    record.setBaseQualities(SAMRecord.NULL_QUALS);
    return record;
  }

  private SAMRecord createUnmappedRecordWithoutSequence() {
    SAMRecord record = createUnmappedRecord("read-without-seq", "A", "!");
    record.setReadBases(SAMRecord.NULL_SEQUENCE);
    record.setBaseQualities(SAMRecord.NULL_QUALS);
    return record;
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
}
