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

import static org.junit.Assert.*;
import static uk.ac.ebi.ena.readtools.refactored.TestFileUtil.createOutputFolder;
import static uk.ac.ebi.ena.readtools.refactored.TestFileUtil.saveRandomized;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import java.util.zip.GZIPOutputStream;

import org.junit.Assert;
import org.junit.Test;

import uk.ac.ebi.ena.readtools.refactored.provider.FastqReadsProvider;
import uk.ac.ebi.ena.readtools.refactored.provider.ReadsProvider;
import uk.ac.ebi.ena.readtools.refactored.validator.FastqReadsValidator;
import uk.ac.ebi.ena.readtools.refactored.validator.ReadsValidationException;

public class FastqReadsValidatorTest {
    @Test
    public void validSingleFile1() throws ReadsValidationException {
        ReadsProvider mrp = new FastqReadsProvider(
                Paths.get("src/test/resources/rawreads/EP0_GTTCCTT_S1.txt.gz").toFile());
        assertTrue(new FastqReadsValidator().validate(mrp));
    }

    @Test
    public void validSingleFile2() throws ReadsValidationException {
        ReadsProvider mrp = new FastqReadsProvider(
                Paths.get("src/test/resources/rawreads/SPOP-87C_plKO2-min.fastq.gz").toFile());
        assertTrue(new FastqReadsValidator().validate(mrp));
    }

    @Test
    public void readNameDuplicate1() throws IOException {
        File output_dir = createOutputFolder();
        Path f1 = saveRandomized( "@NAME1/1\nACGT\n+\n1234\n"
                + "@NAME1/2\nACGT\n+\n1234\n"
                + "@NAME1/1\nACGT\n+\n1234", output_dir.toPath(), true, "fastq-1", "gz" );

        try {
            ReadsProvider mrp = new FastqReadsProvider(f1.toFile());
            new FastqReadsValidator().validate(mrp);
            fail();
        } catch (ReadsValidationException e) {
            assertTrue(e.getMessage().contains("Multiple"));
        }
    }

    @Test
    public void readNameDuplicate2() throws IOException {
        File output_dir = createOutputFolder();
        Path f1 = saveRandomized( "@NAME1/1\nACGT\n+\n1234\n"
                + "@NAME1/1\nACGT\n+\n1234\n"
                + "@NAME2/1\nACGT\n+\n1234", output_dir.toPath(), true, "fastq-1", "gz" );

        try {
            ReadsProvider mrp = new FastqReadsProvider(f1.toFile());
            new FastqReadsValidator().validate(mrp);
            fail();
        } catch (ReadsValidationException e) {
            assertTrue(e.getMessage().contains("Multiple"));
        }
    }

    @Test
    public void readNameDuplicate3() throws IOException {
        File output_dir = createOutputFolder();
        Path f1 = saveRandomized( "@NAME1/1\nACGT\n+\n1234\n"
                + "@NAME1/2\nACGT\n+\n1234\n"
                + "@NAME1/1\nACGT\n+\n1234", output_dir.toPath(), true, "fastq-1", "gz" );

        try {
            ReadsProvider mrp = new FastqReadsProvider(f1.toFile());
            new FastqReadsValidator().validate(mrp);
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

        ReadsProvider mrp = new FastqReadsProvider(f1.toFile());
        assertTrue(new FastqReadsValidator().validate(mrp));
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
            ReadsProvider mrp = new FastqReadsProvider(f1.toFile());
            assertTrue(new FastqReadsValidator().validate(mrp));
        } catch (ReadsValidationException e) {
            assertTrue(e.getMessage().contains("Multiple"));
        }
    }

    @Test
    public void testInvalid() {
        File f1 = Paths.get("src/test/resources/invalid.fastq.gz").toFile();
        try {
            ReadsProvider mrp = new FastqReadsProvider(f1);
            new FastqReadsValidator().validate(mrp);
            fail();
        } catch (ReadsValidationException e) {
            assertEquals(1, e.getReadIndex());
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

        ReadsProvider mrp = new FastqReadsProvider(f1.toFile());
        assertTrue(new FastqReadsValidator().validate(mrp));
    }
}
