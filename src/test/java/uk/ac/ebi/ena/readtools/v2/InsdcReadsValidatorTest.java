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

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import org.junit.Ignore;
import org.junit.Test;
import uk.ac.ebi.ena.readtools.v2.MockReadsProvider.MockRead;
import uk.ac.ebi.ena.readtools.v2.provider.ReadsProvider;
import uk.ac.ebi.ena.readtools.v2.provider.ReadsProviderFactory;
import uk.ac.ebi.ena.readtools.v2.validator.FastqReadsValidator;
import uk.ac.ebi.ena.readtools.v2.validator.InsdcReadsValidator;
import uk.ac.ebi.ena.readtools.v2.validator.ReadsValidationException;

public class InsdcReadsValidatorTest {
  @Test
  public void noReads() {
    try {
      MockReadsProviderFactory factory = new MockReadsProviderFactory();
      new InsdcReadsValidator(READ_COUNT_LIMIT).validate(factory);
      fail();
    } catch (ReadsValidationException e) {
      assertTrue(e.getErrorMessage().contains(ERROR_NO_READS));
    }
  }

  @Test
  public void emptyRead() {
    ReadsProvider mrp;

    try {
      MockReadsProviderFactory factory =
          new MockReadsProviderFactory(new MockRead("r1", null, "1234"));
      new InsdcReadsValidator(READ_COUNT_LIMIT).validate(factory);
      fail();
    } catch (ReadsValidationException e) {
      assertTrue(e.getErrorMessage().contains(ERROR_EMPTY_READ));
    }

    try {
      MockReadsProviderFactory factory =
          new MockReadsProviderFactory(new MockRead("r1", "", "1234"));
      new InsdcReadsValidator(READ_COUNT_LIMIT).validate(factory);
      fail();
    } catch (ReadsValidationException e) {
      assertTrue(e.getErrorMessage().contains(ERROR_EMPTY_READ));
    }

    //whitespace as read
    try {
      MockReadsProviderFactory factory =
              new MockReadsProviderFactory(new MockRead("r1", " ", "1234"));
      new InsdcReadsValidator(READ_COUNT_LIMIT).validate(factory);
      fail();
    } catch (ReadsValidationException e) {
      assertTrue(e.getErrorMessage().contains(ERROR_EMPTY_READ));
    }
  }

  @Test
  public void testZeroReadsSingleFile() throws IOException {
    File output_dir = createOutputFolder();
    Path f1 = saveRandomized("", output_dir.toPath(), true, "fastq-zero-reads-single", "gz");

    try {
      ReadsProviderFactory factory = new ReadsProviderFactory(f1.toFile(), FileFormat.FASTQ);
      new InsdcReadsValidator(READ_COUNT_LIMIT).validate(factory);
      fail();
    } catch (ReadsValidationException e) {
      ;
      // todo: assertEquals(ERROR_NO_READS, e.getMessage());
    }
  }

  @Test
  public void testZeroNucleotidesSingleFile() throws IOException {
    File output_dir = createOutputFolder();
    Path f1 =
        saveRandomized(
            "@CL100031435L2C001R001_63456 :TGGTCCTTCC\n"
                + "\n"
                + "+\n"
                + "FDFFFFFFFG3FFFFFF7?D8FFFFF=>DFG(FBEFFFFDF;FF?FF<8FGGFFFCBB8F0F@FBC?FAGFEE>.FFEFCF:?E(E@;3*(,FD/BFE-\n",
            output_dir.toPath(),
            true,
            "fastq-10",
            "gz");

    try {
      ReadsProviderFactory factory = new ReadsProviderFactory(f1.toFile(), FileFormat.FASTQ);
      new InsdcReadsValidator(READ_COUNT_LIMIT).validate(factory);
      fail();
    } catch (ReadsValidationException e) {
      ;
      // todo: assertTrue(e.getMessage().contains("Missing Sequence Line"));
    }
  }

  @Test
  public void testReadNameLength() throws IOException {
    File output_dir = createOutputFolder();
    Path f1 =
        saveRandomized(
            "@aaaabbbbccccddddaaaabbbbccccddddaaaabbbbccccddddaaaabbbbccccddddaaaabbbbccccddddaaaabbbbccccddddaaaabbbbccccddddaaaabbbbccccddddaaaabbbbccccddddaaaabbbbccccddddaaaabbbbccccddddaaaabbbbccccddddaaaabbbbccccddddaaaabbbbccccddddaaaabbbbccccddddaaaabbbbccccdddd1\n"
                + "AGTC\n"
                + "+\n"
                + "1234\n",
            output_dir.toPath(),
            true,
            "fastq",
            "gz");

    try {
      ReadsProviderFactory factory = new ReadsProviderFactory(f1.toFile(), FileFormat.FASTQ);
      new InsdcReadsValidator(READ_COUNT_LIMIT).validate(factory);
      fail();
    } catch (ReadsValidationException e) {
      assertTrue(e.getMessage().contains("Read name length exceeds 256 characters"));
    }
  }

  @Test
  public void testZeroQualitiesSingleFile() throws IOException {
    File output_dir = createOutputFolder();
    Path f1 =
        saveRandomized(
            "@CL100031435L2C001R001_63456 :TGGTCCTTCC\n"
                + "ATCCTTATAGTGGGCCAAGCTCTCACATGCAAACACGTTGCCTGCTGGATTTTGTTTCAGAGAGGAAATGTTTATGTGAGACAGAAAAAGCCGGGGGCC\n"
                + "+\n",
            output_dir.toPath(),
            true,
            "fastq-10",
            "gz");

    try {
      ReadsProviderFactory factory = new ReadsProviderFactory(f1.toFile(), FileFormat.FASTQ);
      new InsdcReadsValidator(READ_COUNT_LIMIT).validate(factory);
      fail();
    } catch (ReadsValidationException e) {
      assertTrue(e.getMessage().contains("File is too short"));
    }
  }

  @Test
  public void testEmptyReadSingleFile() throws IOException {
    File output_dir = createOutputFolder();
    Path f1 =
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
            output_dir.toPath(),
            true,
            "fastq-10",
            "gz");

    try {
      ReadsProviderFactory factory = new ReadsProviderFactory(f1.toFile(), FileFormat.FASTQ);
      new InsdcReadsValidator(READ_COUNT_LIMIT).validate(factory);
      fail();
    } catch (ReadsValidationException e) {
      ;
      // todo: assertEquals(INVALID_FILE, e.getMessage());
    }
  }

  @Test
  public void notIUPAC() {
    try {
      MockReadsProviderFactory factory =
          new MockReadsProviderFactory(new MockRead("r1", "AFFF", "1234"));
      new InsdcReadsValidator(READ_COUNT_LIMIT).validate(factory);
      fail();
    } catch (ReadsValidationException e) {
      assertTrue(e.getErrorMessage().contains(ERROR_NOT_IUPAC));
    }
  }

  @Test
  public void notAUTCG() {
    try {
      MockReadsProviderFactory factory =
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
      new InsdcReadsValidator(READ_COUNT_LIMIT).validate(factory);
      fail();
    } catch (ReadsValidationException e) {
      assertTrue(e.getErrorMessage().contains(ERROR_NOT_AUTCG));
    }
  }

  @Ignore("rejects too many submissions")
  @Test
  public void lowQuality() throws ReadsValidationException {
    MockReadsProviderFactory factory =
        new MockReadsProviderFactory(new MockRead("r1", "AGTC", ">>>>"));
    new InsdcReadsValidator(READ_COUNT_LIMIT).validate(factory);

    try {
      factory = new MockReadsProviderFactory(new MockRead("r1", "AGTC", "@@@@"));
      new InsdcReadsValidator(READ_COUNT_LIMIT).validate(factory);
      fail();
    } catch (ReadsValidationException e) {
      assertEquals(ERROR_QUALITY, e.getErrorMessage());
    }
  }

  @Test
  public void validFastq() throws IOException, ReadsValidationException {
    File output_dir = createOutputFolder();
    Path f1 =
        saveRandomized(
            "@CL100031435L2C001R001_63456 :TGGTCCTTCC\n"
                + "ATCCTTATAGTGGGCCAAGCTCTCACATGCAAACACGTTGCCTGCTGGATTTTGTTTCAGAGAGGAAATGTTTATGTGAGACAGAAAAAGCCGGGGGCC\n"
                + "+\n"
                + "FDFFFFFFFG?FFFFFF??D?FFFFF??DFG?FBEFFFFDF?FF?FF??FGGFFFCBB?F?F@FBC?FAGFEE??FFEFCF??E?E@?????FD?BFE?\n"
                + "@CL100031435L2C001R001_63487 :TGGTCCTTCC\n"
                + "CCGCTCTGCCCACCACAATCCGGCCCCTGTGTACGGCAACACAGGGCCCAGGCAAGCAGATCCTTCCTGCTGGGAGCTCCAGCTTGTAGAATTTCACCC\n"
                + "+\n"
                + "FFFFFFFFFEFFFFFFD?FEFFEEEFFEBCE?EEFFFECD?DAFFCDFFBFFDB@?FCFBEEFF?F?CBFFF?E?F?BBF??@B?B??E?F?A?EFFCA\n",
            output_dir.toPath(),
            true,
            "fastq-10",
            "gz");

    ReadsProviderFactory factory = new ReadsProviderFactory(f1.toFile(), FileFormat.FASTQ, true);
    new InsdcReadsValidator(READ_COUNT_LIMIT).validate(factory);
    new FastqReadsValidator(READ_COUNT_LIMIT).validate(factory);
  }
}
