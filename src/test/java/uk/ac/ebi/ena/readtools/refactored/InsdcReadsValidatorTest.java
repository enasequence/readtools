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
package uk.ac.ebi.ena.readtools.refactored;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;
import static uk.ac.ebi.ena.readtools.refactored.validator.InsdcReadsValidator.*;

import org.junit.Test;

import uk.ac.ebi.ena.readtools.refactored.MockReadsProvider.MockRead;
import uk.ac.ebi.ena.readtools.refactored.provider.ReadsProvider;
import uk.ac.ebi.ena.readtools.refactored.validator.InsdcReadsValidator;
import uk.ac.ebi.ena.readtools.refactored.validator.ReadsValidationException;

public class InsdcReadsValidatorTest {
    @Test
    public void noReads() {
        try {
            ReadsProvider mrp = new MockReadsProvider();
            new InsdcReadsValidator().validate(mrp);
            fail();
        } catch (ReadsValidationException e) {
            assertEquals(ERROR_NO_READS, e.getErrorMessage());
        }
    }

    @Test
    public void emptyRead() {
        ReadsProvider mrp;

        try {
            mrp = new MockReadsProvider(
                    new MockRead("r1", null, "1234"));
            new InsdcReadsValidator().validate(mrp);
            fail();
        } catch (ReadsValidationException e) {
            assertEquals(ERROR_EMPTY_READ, e.getErrorMessage());
        }

        try {
            mrp = new MockReadsProvider(
                    new MockRead("r1", "", "1234"));
            new InsdcReadsValidator().validate(mrp);
            fail();
        } catch (ReadsValidationException e) {
            assertEquals(ERROR_EMPTY_READ, e.getErrorMessage());
        }
    }
//
//    @Test public void
//    testZeroReadsSingleFile() throws IOException {
//        File output_dir = createOutputFolder();
//        Path f1 = saveRandomized(
//                "",
//                output_dir.toPath(), true, "fastq-zero-reads-single", "gz");
//
//        ReadsProvider mrp = new FastqReadsProvider(f1.toFile());
//        try {
//            new FastqReadsValidator().validate(mrp);
//            fail();
//        } catch (ReadsValidationException e) {
//            assertEquals("", e.getMessage());
//        }
//    }
//
//    @Test public void
//    testZeroNucleotidesSingleFile() throws IOException {
//        File output_dir = createOutputFolder();
//        Path f1 = saveRandomized(
//                "@CL100031435L2C001R001_63456 :TGGTCCTTCC\n"
//                        + "\n"
//                        + "+\n"
//                        + "FDFFFFFFFG3FFFFFF7?D8FFFFF=>DFG(FBEFFFFDF;FF?FF<8FGGFFFCBB8F0F@FBC?FAGFEE>.FFEFCF:?E(E@;3*(,FD/BFE-\n",
//                output_dir.toPath(), true, "fastq-10", "gz");
//
//        ReadsProvider mrp = new FastqReadsProvider(f1.toFile());
//        try {
//            new FastqReadsValidator().validate(mrp);
//            fail();
//        } catch (ReadsValidationException e) {
//            assertEquals("", e.getMessage());
//        }
//    }
//
//    @Test public void
//    testZeroQualitiesSingleFile() throws IOException {
//        File output_dir = createOutputFolder();
//        Path f1 = saveRandomized(
//                "@CL100031435L2C001R001_63456 :TGGTCCTTCC\n"
//                        + "ATCCTTATAGTGGGCCAAGCTCTCACATGCAAACACGTTGCCTGCTGGATTTTGTTTCAGAGAGGAAATGTTTATGTGAGACAGAAAAAGCCGGGGGCC\n"
//                        + "+\n",
//                output_dir.toPath(), true, "fastq-10", "gz");
//
//        ReadsProvider mrp = new FastqReadsProvider(f1.toFile());
//        try {
//            new FastqReadsValidator().validate(mrp);
//            fail();
//        } catch (ReadsValidationException e) {
//            assertEquals("", e.getMessage());
//        }
//    }
//
//    @Test public void
//    testEmptyReadSingleFile() throws IOException {
//        File output_dir = createOutputFolder();
//        Path f1 = saveRandomized(
//                "@CL100031435L2C001R001_63456 :TGGTCCTTCC\n"
//                        + "ATCCTTATAGTGGGCCAAGCTCTCACATGCAAACACGTTGCCTGCTGGATTTTGTTTCAGAGAGGAAATGTTTATGTGAGACAGAAAAAGCCGGGGGCC\n"
//                        + "+\n"
//                        + "FDFFFFFFFG3FFFFFF7?D8FFFFF=>DFG(FBEFFFFDF;FF?FF<8FGGFFFCBB8F0F@FBC?FAGFEE>.FFEFCF:?E(E@;3*(,FD/BFE-\n"
//                        + "@CL100031435L2C001R001_63472 :TGGTCCTTCA\n"
//                        + "\n"
//                        + "+\n"
//                        + "\n"
//                        + "@CL100031435L2C001R001_63487 :TGGTCCTTCC\n"
//                        + "CCGCTCTGCCCACCACAATCCGGCCCCTGTGTACGGCAACACAGGGCCCAGGCAAGCAGATCCTTCCTGCTGGGAGCTCCAGCTTGTAGAATTTCACCC\n"
//                        + "+\n"
//                        + "FFFFFFFFFEFFFFFFD5FEFFEEEFFEBCE?EEFFFECD:DAFFCDFFBFFDB@>FCFBEEFF>F-CBFFF&E9F=BBF:4@B?B45E2F+A-EFFCA\n",
//                output_dir.toPath(), true, "fastq-10", "gz");
//
//        ReadsProvider mrp = new FastqReadsProvider(f1.toFile());
//        try {
//            new FastqReadsValidator().validate(mrp);
//            fail();
//        } catch (ReadsValidationException e) {
//            assertEquals("", e.getMessage());
//        }
//    }

    @Test
    public void notIupac() {
        try {
            ReadsProvider mrp = new MockReadsProvider(
                    new MockRead("r1", "AFFF", "1234"));
            new InsdcReadsValidator().validate(mrp);
            fail();
        } catch (ReadsValidationException e) {
            assertEquals(ERROR_NOT_IUPAC, e.getErrorMessage());
        }
    }

    @Test
    public void lowQuality() throws ReadsValidationException {
        ReadsProvider mrp = new MockReadsProvider(
                new MockRead("r1", "AGTC", ">>>>"));
        new InsdcReadsValidator().validate(mrp);

        try {
            mrp = new MockReadsProvider(
                    new MockRead("r1", "AGTC", "@@@@"));
            new InsdcReadsValidator().validate(mrp);
            fail();
        } catch (ReadsValidationException e) {
            assertEquals(ERROR_QUALITY, e.getErrorMessage());
        }
    }
}
