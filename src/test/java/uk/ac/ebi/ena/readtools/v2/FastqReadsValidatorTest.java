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

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Collections;

import org.junit.Test;

import uk.ac.ebi.ena.readtools.v2.provider.FastqReadsProvider;
import uk.ac.ebi.ena.readtools.v2.provider.ReadsProvider;
import uk.ac.ebi.ena.readtools.v2.provider.ReadsProviderFactory;
import uk.ac.ebi.ena.readtools.v2.validator.FastqReadsValidator;
import uk.ac.ebi.ena.readtools.v2.validator.ReadsValidationException;
import uk.ac.ebi.ena.readtools.v2.validator.ValidatorWrapper;

public class FastqReadsValidatorTest {
    @Test
    public void validSingleFile1() throws ReadsValidationException {
        ReadsProviderFactory factory = new ReadsProviderFactory(
                Paths.get("src/test/resources/rawreads/EP0_GTTCCTT_S1.txt.gz").toFile(), FileFormat.FASTQ);
        assertTrue(new FastqReadsValidator(READ_COUNT_LIMIT).validate(factory));
    }

    @Test
    public void validSingleFile2() throws ReadsValidationException {
        ReadsProviderFactory factory = new ReadsProviderFactory(
                Paths.get("src/test/resources/rawreads/SPOP-87C_plKO2-min.fastq.gz").toFile(), FileFormat.FASTQ);
        assertTrue(new FastqReadsValidator(READ_COUNT_LIMIT).validate(factory));
    }

    @Test
    public void readNameDuplicate1() throws IOException {
        File output_dir = createOutputFolder();
        Path f1 = saveRandomized( "@NAME1/1\nACGT\n+\n1234\n"
                + "@NAME1/2\nACGT\n+\n1234\n"
                + "@NAME1/1\nACGT\n+\n1234", output_dir.toPath(), true, "fastq-1", "gz" );

        try {
            ReadsProviderFactory factory = new ReadsProviderFactory(f1.toFile(), FileFormat.FASTQ);
            new FastqReadsValidator(READ_COUNT_LIMIT).validate(factory);
            fail();
        } catch (ReadsValidationException e) {
            assertTrue(e.getMessage().contains("Multiple"));
        }
    }

    @Test
    public void readNameDuplicate2() throws IOException {
        File output_dir = createOutputFolder();
        Path f1 = saveRandomized(
                "@NAME1/1\nACGT\n+\n1234\n" +
                        "@NAME1/1\nACGT\n+\n1234\n" +
                        "@NAME2/1\nACGT\n+\n1234",
                output_dir.toPath(), true, "fastq-1", "gz" );

        try {
            ValidatorWrapper vw = new ValidatorWrapper(
                    Collections.singletonList(f1.toFile()), FileFormat.FASTQ, READ_COUNT_LIMIT);
            vw.run();
//            ReadsProvider mrp = new FastqReadsProvider(f1.toFile());
//            new FastqReadsValidator(READ_COUNT_LIMIT).validate(mrp);
            fail();
        } catch (ReadsValidationException e) {
            assertTrue(e.getMessage().contains("Multiple"));
        }
    }

    @Test
    public void readNameDuplicate3() throws IOException {
        File output_dir = createOutputFolder();
        Path f1 = saveRandomized(
                "@NAME1/1\nACGT\n+\n1234\n"
                + "@NAME1/2\nACGT\n+\n1234\n"
                + "@NAME1/1\nACGT\n+\n1234\n"
                + "@NAME1/1\nACGT\n+\n1234",
                output_dir.toPath(), true, "fastq-1", "gz" );

        try {
            ReadsProviderFactory factory = new ReadsProviderFactory(f1.toFile(), FileFormat.FASTQ);
            new FastqReadsValidator(READ_COUNT_LIMIT).validate(factory);
            fail();
        } catch (ReadsValidationException e) {
            assertTrue(e.getMessage().contains("Multiple"));
        }
    }

    @Test
    public void validPairedReads1() throws IOException, ReadsValidationException {
        File output_dir = createOutputFolder();
        Path f1 = saveRandomized(
                "@m150522_164945_42237_c100824932550000001823179511031571_s1_p0/134/0_15617\n"
                        + "ACCAACAAGAGA\n"
                        + "+\n"
                        + "\"#'#&''$,-,)\n"

                        + "@m150522_164945_42237_c100824932550000001823179511031571_s1_p0/135/0_15617\n"
                        + "ACCAACAAGAGA\n"
                        + "+\n"
                        + "\"#'#&''$,-,)\n"

                        + "@m150522_164945_42237_c100824932550000001823179511031571_s1_p0/134/15661_24597\n"
                        + "ACCACCCTTAT\n"
                        + "+\n"
                        + "/&\"-('--.#/\n"

                        + "@m150522_164945_42237_c100824932550000001823179511031571_s1_p0/135/15661_24597\n"
                        + "ACCACCCTTAT\n"
                        + "+\n"
                        + "/&\"-('--.#/\n",
                output_dir.toPath(), true, "fastq-9", "gz" );

        ReadsProviderFactory factory = new ReadsProviderFactory(f1.toFile(), FileFormat.FASTQ);
        assertTrue(new FastqReadsValidator(READ_COUNT_LIMIT).validate(factory));
    }

    @Test
    public void validPairedReads2() throws IOException, ReadsValidationException {
        File output_dir = createOutputFolder();
        Path f1 = saveRandomized(
                "@CL100031435L2C001R001_63456 :TGGTCCTTCC\n"
                        + "ATCCTTATAGTGGGCCAAGCTCTCACATGCAAACACGTTGCCTGCTGGATTTTGTTTCAGAGAGGAAATGTTTATGTGAGACAGAAAAAGCCGGGGGCC\n"
                        + "+\n"
                        + "FDFFFFFFFG3FFFFFF7?D8FFFFF=>DFG(FBEFFFFDF;FF?FF<8FGGFFFCBB8F0F@FBC?FAGFEE>.FFEFCF:?E(E@;3*(,FD/BFE-\n"
                        + "@CL100031435L2C001R001_63472 :TGGTCCTTCA\n"
                        + "AGCACTATGGAAAATGGAACCTCTTTGGGTCCTGTAATCAC\n"
                        + "+\n"
                        + "FFFFFFEGFFFGEFGFFFGFGGFEFDFFFFGGEFEFBG;FG\n"
                        + "@CL100031435L2C001R001_63487 :TGGTCCTTCC\n"
                        + "CCGCTCTGCCCACCACAATCCGGCCCCTGTGTACGGCAACACAGGGCCCAGGCAAGCAGATCCTTCCTGCTGGGAGCTCCAGCTTGTAGAATTTCACCC\n"
                        + "+\n"
                        + "FFFFFFFFFEFFFFFFD5FEFFEEEFFEBCE?EEFFFECD:DAFFCDFFBFFDB@>FCFBEEFF>F-CBFFF&E9F=BBF:4@B?B45E2F+A-EFFCA\n"
                        + "@CL100031435L2C001R001_63499 :TGGTCCTTGC\n"
                        + "CTCTTTTGATGCCTTCTGTATTGAGACACTTTCAGCAGAGTCCCAAGGCCAGGTGAAGAGGAAGCATGGCGGGGGCAGGTCGGGGGCCAACGGGTCTTAG\n"
                        + "+\n"
                        + "FEFF@E9F=FFF7FDD@FDFFFF4F;EFFFFAA8FEADDFFEFFE5F?FEBDEFFACD;BFFBFD2F>EDD:7FDDF'FDFF)+7FF1F1<4FF:@A@/F\n"
                        + "@CL100031435L2C001R001_63541 :TGGTCCTTGC\n"
                        + "CATGGGGGCATCTTCCTGCTCCAGGCACAAAGCTCTGGCTGTCACAACCCAG\n"
                        + "+\n"
                        + "FFFGFF>GFFGFFBFAGGFEFFDEFGFFEEAGGFFFBB>GGBFF6FFFGF@F\n"
                        + "@CL100031435L2C001R001_63542 :TGGTCCTTGC\n"
                        + "CCTCATTCTAAAATGCAGGTGCTCTGCTGACAGCAAAATTCTGTGTTTGAGCTGTGCAGTCTTAAACCAAGACGTGGGAGGCCAGGCGGTTTTAGAAAG\n"
                        + "+\n"
                        + "AAE>AFE4EFFA3FBAEDDFBEFDDA?F=;E?6?>FC4BF@EFC5A9FE=ECFCCD<CEEBF(EACE1F<5D(0(+DA;C9@1F&3,F76&B?DF1<30\n"
                        + "@CL100031435L2C001R001_63548 :TGGTCCTTGC\n"
                        + "GTGGTGGTGGTGTTTTTTTGGTTTGGGGTTTGGGTTTTTTTCAGGGAA\n"
                        + "+\n"
                        + "GFGFFFFFF=DGFFGGEFFFFFBGFF>FFFGFFFEGFGGEEGGGGFGF\n"
                        + "@CL100031435L2C001R001_63553 :TGGTCCTTCC\n"
                        + "GCACTTCTCGAGCTTCACATTCTAATGAGAACAATTTCCTTGGATTCATTGGTGTTGCCATTTTTTTGTTGACTCATTCAAAAAACAAATTAGCTGAGGT\n"
                        + "+\n"
                        + "FFFFFFFFFFDFEGFFGEGG=EEDFGFEFEFEEGEFFFFEFAD:BEDGBGFDF,FFAF7<FF;(GEGFFD8?FBF?/1-GGD8'FF'ECGA:AF<FF&AD\n"
                        + "@CL100031435L2C001R001_63565 :TGGTCCTTGC\n"
                        + "AAAGCTCTAGGGGAGGCCGACCTCTCAGCTTTTGGAGTCGGTGATACAGAGGA\n"
                        + "+\n"
                        + "EFFEDBCFFDAC7C@DBEF6DDFDFE;EEFFF:DA=DFDE8EC;DF@E0=AEF\n"
                        + "@CL100031435L2C001R001_63596 :TGGTCCTTGC\n"
                        + "CCTTCTTTTTGATCCTACAAGTGAAGCCATGTGGAAATGATTGAGAGTGACGTCACGAGTTCAGTTGTGGGGGGCGGCGGCCGGGTTGTCCGGTCGGGAC\n"
                        + "+\n"
                        + "CBEC>=EFCED@?CDDF<EFA>4CDD8E8AAC>:4FDD>=CDD4F7CEF;CAB=FB>:CAA?F,EFE?;&@D7>5:?/==(C>.<&26'-=?8%-+E42F",
                output_dir.toPath(), true, "fastq-10", "gz");

        try {
            ReadsProviderFactory factory = new ReadsProviderFactory(f1.toFile(), FileFormat.FASTQ);
            assertTrue(new FastqReadsValidator(READ_COUNT_LIMIT).validate(factory));
        } catch (ReadsValidationException e) {
            assertTrue(e.getMessage().contains("Multiple"));
        }
    }

    @Test
    public void testInvalid() {
        File f1 = Paths.get("src/test/resources/invalid.fastq.gz").toFile();
        try {
            ReadsProviderFactory factory = new ReadsProviderFactory(f1, FileFormat.FASTQ);
            new FastqReadsValidator(READ_COUNT_LIMIT).validate(factory);
            fail();
        } catch (ReadsValidationException e) {
            assertEquals(0, e.getReadIndex());
        }
    }

    @Test
    public void testSinglePairedWithUnpairedFastq() throws IOException, ReadsValidationException {
        File output_dir = createOutputFolder();
        Path f1 = saveRandomized(
                "@paired/1\n" +
                        "ACGT\n" +
                        "+\n" +
                        "1234\n" +
                        "@paired/2\n" +
                        "ACGT\n" +
                        "+\n" +
                        "1234\n" +
                        "@unpaired\n" +
                        "ACGT\n" +
                        "+\n" +
                        "1234", output_dir.toPath(), true, "fastq-1", "gz");

        ReadsProviderFactory factory = new ReadsProviderFactory(f1.toFile(), FileFormat.FASTQ);
        assertTrue(new FastqReadsValidator(READ_COUNT_LIMIT).validate(factory));
    }
}
