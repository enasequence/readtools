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
import java.util.Arrays;
import java.util.Collections;

import org.junit.Ignore;
import org.junit.Test;

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
    public void validSingleFile3() throws IOException, ReadsValidationException {
        File output_dir = createOutputFolder();
        Path f1 = saveRandomized(
                "@Chr3 (41795662..41800467) consensus\n" +
                        "ATTCCAAACGCAGTTTTATTTCATACAAGAATTATACAATATAAAATCATTCCCTTTTTCCTCTCTTTGTGAATAGAAACTACATATATCAGAAAGCCACTTTCCGGACACAGTCACAACCCAAAAAACCCAAATGCAATAGTTATGTACCCAAGCTGAGGACTCAAATCTTTTTAAGGAATTCCATGATGTAAGTGTTGAAAGCAGAACACACTACTTTGAACCCTTGAAATCCTGCCCCTTTAGCCAAGTCCTCGAACTCCTTCTCGGTCCTCTCTTTCCCACCCGGATTATGAGCCAACATGATAACATCAATATGGATGACTCCCTTAGCAGCAAGGCTACTGTCCGGGGCTACAGGAAGAATACATTCAGCAACAATCACCTTCCCATTTTCTGGAAGTGCATCATAGCAGTTCTTCAAAAATTTTGAGCAGTGCTCATCACTCCAATCATGGCATATCCACTGTACATAGGAGCACACACATGTCAGTAGCATAATCAAAGATCAAATGTCTTTGACTAATCAATTTTTCTTTAAATTTCATTTTCACCTTATGATTCATGTGCTAGTTAATAGGTTAAAGTGTCTTTCCACAAAAAGTATACTCTTCTAGTTACTATTATAAGCACAAACACTATCTCTCCACTTAGCTTCCCCATCCCTTAAAATGGATTTTTTAACAATACACTGATATCACTTTAGTGCTGCATTAAAATTGGTTTCCCTTTAAGTTTTTTTAAAAATGAAACACACCCTATGGTTGAATTAGAGAGTAATCCTGACTCCTGAGTGCTAGCAGTGCCTTCCAAAGTGTCAACTATCAATATCATCATTGTAGGTTTATATAACTGGACCCAACTAAAAGTCAAGCAAGCGTGTCGGTAAGAAAGTGAAGGGATTTCAGATAAGTTCAAACCCACCTTCATGAAAATGGCATCTCCTTTAGGAACACTAACAAACATGTCTCCCCCAACATGCTCCACACCTGCAGCCAACACCCATTTCTAGCCCATGAGATATTTAACCCAAAAACATAAAAAAGGAACAACAATATTATTATGTAGCCCACATGGTTTAGTTTGGGAGTCAAATAATACTCTTTAAACAAAAAAAGACTCCAAAAGAGACCACTGGTCTATATCTATGAAAATATTTAGTGAACCCAGTGTACTACACATTTTTCTCACATAAAAGTAGACTTTATATGTGAATCCCGCATTAGGACCCATAAGTGTAGTATACTGGGTGTACAAAATAATTTTTCACATCCAACAACCTCTAAAAGAGTCCAACTTATTAAAGAGAAAGGTCCAAATAGGCAAGCACTTAAGTCTACTTCACCATTATAAAATTAAGCATGCTATTAATCACATGACTTAAGTAAACACAAAATAATGCCAATACCAGGATAAGATGGTGCATCCTCAATAACATGGGGCAAATCAAAATTGATGCCCTTTATAGAGGGGTATTTCGAGACAATCGTATTAATGATAGCTCCAGTCCCACCACCAACATCAACCACAGAAGTGAGACCCTCGAAGCCTTTGTAGGTCTCCAATATTTTCTTCATGGTAATGGTAGAGTGATCAGACATTCCCTTGTTGAAGACCTTGTTGAATCTAAGATCTTTCCCATGATACTCAAAAGAGGTCATTCCATGGGCCTTGTTAAATGGAATGCCACCTTCAAGAACTGCCTCTTTCAAGTAGTACCTAAAGAAAACATTATTAGATAAAAAATTATTAGCAACCAAAAATACAAACCTAAGAAATCAAGTGTGAAATCGTGGAATTTTACAAATTGTGATGAGTTCTTTTCTCATCATTGTATAACCTTTCTCGGGGTGGTTGTTACCAGCACCCCGAGAAAGGTTAGTCCAGCCAAAAATAATTATTAGCAACCAAAAATACAAACCAAATAATTCTAGATTTGTTTTTCAATTTTTCTTTCCTCTAATCCTTCTTCCAGAAAACATATTTTATGTATCGGATATATATTCATGAGATAATTATTAATGAGATGAATATATATATATATTATAATTGCCCATAGGAATCATGTGGCTTCCAATTGTTGAGATATGATAATTAGATGCATAGAATATAGAAGTTATGCACAGCCTTAGTAAAGAAGTTAGGTTGAGTTATCGTCTAAAAACAAAAAAGTAAAGCAAAAAGGAGATGATAGGTGAAGTGGTGGTTGATCATTGATGGTTCTCGAGTCCTAAGTGGACGGTACCCCATGCCTAGCCCATCAGCAAGTTGTGACTATCTATCGCACAAATTACCACATAAACGTTTTAGGGGTCCCACTAAGGAATTTTTATTTTTTATTTTATTTTAATTTTCTTTGCATAATGCGAGTAAAGAAAGAGTAAGAAACAAAGTATGTCTTCAAATTACACCAAACATTCAATTAATAACAAGCATGATAGCCATTGGGTATTCGGCCCAGTAGTAATCTCTCTTACCAAAAAGCTTTACTCCTCTATGAGAGTGTGAGTTCAATCCCCGAAGTCTTCAAATATATGTATATGTGAATTTAGATTAGAATTTCTATCATATATTAAAAATAAAAAAATAAATAACAAGCATGATAATCACCGATGCGGTTGTTGAAGCAACTATTAGGCACTCAAAATGCTTGCTTAGTAAGGCAATAGTACACTCCAATTGGCTGTAAATTAGATTCTCCAATTATATTTTAAATAAATAAATAAAAATTGCTCTCACAATTCTCCTTGATAATTAATATTTTCTGTATCTCTTTTTTTATTGTATTCTAAAATATAATTCTATTTGAAAACTCAATAAACATAATATTAATAACTTATCTATTTTATTTGAAAAATGAATAACCTTTTAAATTAGTTTTATTTGAATTTAAATTAAACAAGATAATTTTAAAAATTTATATTTTTTCAATTGAAAGACAATAAACTTAAAAGCCCAATCTCATTATTTTACTGATTTTTTATAGCCTTGCGCAGAAATTCTCGAAAAGTTATTAGTTCCATCACTTTTATGAAGGTCGTATTATATATTTATATGTAATGATTACATACAACATCTTTCTCACATGGATTGAATCCCACTTTTGTGAGGTTGATGAGATGTTACATGCAACATTTTACATAAGATGTTAAATCCAGTTTACGAGGTGATATCTCATGTTTGAGCTACATGCAGTGCCACTCCACAACAAGCGAGTTCCACGAAACATAAGATAATTTCTCCCAATTGACAGGTGGTATATTATGTAACCGTTACATATAACATGTGTCTCACATGGGGTGAGTCTCATGCTTGTGAGATCCACATCCATATAAGAGAGATGTTAAATATATGGACATAAAATACTTTAAAGTTTACAGAGTCTACAGTGCTCATTCGTGGAACCCACAAATACACAAAATTCTTTTCCGATAGACTCTGCGCATTTGTTTTGATAACATTAGCGTGGGCCTCACTCGTGGATCCCATGAGGCGGCCCAGCGTAATAATCACTCGAAGTCCAGGCAGTAAACGTCTGTACTTGTGATAAATCATTAATGGGATGAAGAAATACTCAAACCATACTCTACATAATTCAAACTTCTTCACAGCCTAGAAAGTACTAAGTAGTAAGAACAAGGCAATAGTACTTACAAAAAGCAATTAATTTTTTTTGTGTGAAAATTTATGAGCTTTTTGTTTAATCCTATGATTTCATTTCATAAGTATAAAATTAATAAATTAGTTTGGGCCGCAATAATAGGAGAAGGCTTCTCTGTTTTAATTAAAGAGAATGACAGTGGTCATAAGAAAATAAAACTTAGAAAAGAAGTGTTGCATGTAATTATTATTACCTATTAGCACAATTAGTAATAGATTTTGTTAACGGGTGTCCTTAGGATAATGGTTAAGGAGTATTAATTATTTACACTTGATGGAAAATAACAAATTTTTAATATTAAAATAATAAATAAAAAGACTGAAATATTTTAAAAATTGTTTTTATATTATTTTAATATTAAAATTTCACTACTTTCTATAAAAATATTAAATATTTTCTTAACCATTGTACTAAGGCCACTCCTTAAAAAGACCTATTAATAATTAACTTCTACTCCAAAAGTTTATAAAACACAAAATTAACTACAAACTAATTTACAAGTTTGTAATTACTACGGATCTTATACCATCCACGTGTATTTTTTAGTCAATCTTATGAGTTGCATGTAACAAAATTACTACAAGTAGTTAAACTTTGTTGTAAAATTAATTAATTAATTAATTAAACCAATAGGAAAGGAAAGAAGAGGAGAATAGATGATTACCAGCTCTCCATGAGGACCTTATCCTGATTCATGAGACAAAGAGGGGCAATAGACACACCATCTTCATTCTTGGTCAAGAATTTACACACAGGGCCTGCACCATAAAGCCTCTCAATCCTGCCATTAGGGAGTGTGCGTAGAGAGTATGTGAGAATAGAGTAGCTAGTCAAGAGACGCAATATACGGTCCACAATTACTGGCGCATCTGGGTTGCTTGTAGGGAGCTGGGAAGCGATTTCTGAAGGTGAAAGGTACGCGCCTGGACCAGCTTTAGCTATGATTTCCAAGAGGTCAAGCTCTATGGCTGATTTGAGAACCATGGGGAGTACTGAGGCACTGGCTAATTGCATGGCAAAGAGGTTTGCTTCTTCATCTGAGACTTGGGTTGGAGTCATCTGGGTTTCACCAGTTGAGCCCATTTTTGTGGTTTGGTTTTGTGGGTTTCTAGTGTGAAAGAGAGAGAAAAGATTAGGAAGATGGG\n" +
                        "+\n" +
                        "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa!aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa5aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa#aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\"aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\n" +
                        "\n" +
                        "",
                output_dir.toPath(), true, "fastq-1", "gz");

        ValidatorWrapper vw = new ValidatorWrapper(
                Arrays.asList(f1.toFile()), FileFormat.FASTQ, READ_COUNT_LIMIT);
        vw.run();
    }

    @Test
    public void readNameDuplicate1() throws IOException {
        File output_dir = createOutputFolder();
        Path f1 = saveRandomized("@NAME1/1\nACGT\n+\n1234\n"
                + "@NAME1/2\nACGT\n+\n1234\n"
                + "@NAME1/1\nACGT\n+\n1234", output_dir.toPath(), true, "fastq-1", "gz");

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
                output_dir.toPath(), true, "fastq-1", "gz");

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
                output_dir.toPath(), true, "fastq-1", "gz");

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
                output_dir.toPath(), true, "fastq-9", "gz");

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
    public void testInvalid2() throws IOException {
        File output_dir = createOutputFolder();
        Path f1 = saveRandomized(
                "@Chr4 (34221604..34224117) consensus\n" +
                        "ACAGCACACAGCACAGTGTTTCTTTCATCACAATGGACCGCTCTTTTCTTTTTAGCCTCTTGATCTTATTCTCTGTGGTCGCGTCCACCGTATTCGCTGATGAGTTAATCAACGGCGATGATGATCCTTTGATCCAACAGGTGGTTTCGGGCGGTGATGATGATCTTCTTCTCCACGCGGAGCATCATTTCTCGAACTTCAAGGCGAAATTCGGAAAGAGTTATGGTAGCAAAGAGGAGCATGATTATAGGCTCGGTGTGTTTAAGTCGAATCTTCGTCGAGCCAAGATGCATCAGAAGCTGGACCCTAACGCAGTCCACGGCGTCACGAAATTCTCCGATCTCACTCCCGCTGAGTTCCGGCGAAACTTTCTTGGCTTGAAGAAGGTTCGGCTCCCGAACGATGCTCAGAAGGCTCCTATTCTTCCCACCAACGATCTTCCCACTGACTTTGACTGGCGCGATCACGGCGCCGTTACCGGCGTCAAAGACCAGGTTTTAAATTTTCAATATAATCAATAATTATATTATTATTATTGTTTTGTTAATTTATAAAGTTTAATCTAAAAATTGAGATTTCATATATTTATATAAATTGTTGGGAATTAGGGTGCGTGTGGATCGTGCTGGTCGTTTAGTACAACGGGAGCTTTGGAAGGAGCTCATTATTTGGCGACTGGGGAGCTAGTGAGCCTCAGTGAACAACAGCTTGTGGATTGTGACCATGAGGTTTGAGTTTGACTCGTTCTTCCTTTGTG   TTGTTTGTTTGGTTTTTTGGTGTGGAGGTGAAATTGTCGTTATGGTTTTGTTGGTCTTTTATTGTTGATAGACTAATTATGTGGTTGTTTTTGTTAGAACAGTAGATAATTGTCTAAATATACACTTTCTTATTAGTTTAAGCTTTTAGGATAGAGCAGATGGTTTAAAAGTTAAAACAAATTGCTTTCATTAGCTTTTGTGAATGAATTTTATTTGAGATTCTTCTACTTTACCTAGTACGTGGCAAAGGTGCCAATTGCCTAAGATACCCAGAAGTGCTCACCTTGTTAATCCCTAGAACTGGAAGTCAAAACTCAAGAGTAATGCTACCTACATTAACATAATTCACAACAAATTTCACAACTGCTGAGATGACAGTTTTTTATTAGTTCTCATCTAGACCTACCACTGACATCACTTTATTATTTATCGTTTATAAGAGCAAATTGTGAAGAAGAA    AAAAAATTCTGAAATTTTGTGTGCCTCTAAACTTATTGAAACAAAATTGAAAGCTTGGGCATGTAGTTGTTCACTATATCAGTATTAGAGATGACTTGCATATAGTTATTTTGTATCTAACCCATGATCATGTGTCGCTATTCTGAAGTTTGACAACCTCGCTAGCGGGAGAGAAAAACCAAAATATGACTATTTCTAAGCATCTTTCTTGCCTCACATTGATTGCTTTGGTCTTGTCTTAGATCACCCCTATGATTCTTTTTGCAGTGTGATCCCGAGGAATATGGTGCATGCGACTCGGGGTGTAATGGTGGGCTGATGACCTCTGCCTTTGAGTACACCCTCAAGGCTGGTGGACTCGAGAGAGAGAAGGATTATCCTTACACCGGGACTGATCGTGGGACCTGCAAATTTGACAAGAGCAAGATTGTTGCTTCTGTATCTAACTTCAGTGTCGTTTCGATTGATGAAGATCAAATTGCCGCGAAATTTGGTGAAGAATGGCCCTCTTGCAAGTAATTACTAATTAACCTTCTTCGAAACCTAGTAAATTGTTTTCATTTCCATCTTTTATTTTTTTTCTTGTTTGGTTGTGTGCTGATGTTGAAATTAGGATAGGCTATGTAACTGTTTTCTTTTTTATTTGAAAAATGGGGTCAAGAATGATAGAGTAAATCGAGTGATGGATATGATATAAATAATGATTCTTTTGCTCTGCTTTGTGCAGTTGGCATCAATGCAGTTTTCATGCAGACCTACATGAAGGGAGTTTCATGCCCATACATCTGCGGGAGGCATTTGGATCATGGTGTGCTTCTGGTGGGTTTTGGATCTGCTGGTTATGCTCCTATCCGGCTCAAGGAGAAGCCTTTCTGGATCATAAAGAACTCCTGGGGAGAAAACTGGGGAGAGAATGGATATTACAAGATCTGTAGGGGCCGTAATGTATGTGGAGTGGATTCCATGGTCTCAACTGTGGCTGCTATAAATGATTAGATTAAGAAGTTTCATGGCAGTGCTGTGTTAGGCAAGTGGAGTCGTTGTAAATATTTAAGTGATAATATGTATTAAGATGACTTCCTAAACTTGTATAAGCTCTTGATGCTTCCATCAAATTACCCAAAAACATGTGTGTCTTGCTCTTGATGGTTCTTTGAAACTTTGTGTCAAGCATTGTCCTGTGTTATTATTCATATTTATCATATGGGTAAGTGAAGCTTGATGCCATGCAGTGGGATGGTTAATTATGTTTTGCTGCGTCCTGATTTGGATTGAAAAATAATGAGTTGATGTACTA\n" +
                        "+\n" +
                        "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa!aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\n",
                output_dir.toPath(), true, "fastq-1", "gz");

        try {
            ValidatorWrapper vw = new ValidatorWrapper(
                    Arrays.asList(f1.toFile()), FileFormat.FASTQ, READ_COUNT_LIMIT);
            vw.run();
            fail();
        } catch (ReadsValidationException e) {
            e.printStackTrace();
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

    @Test
    public void
    testValidPair() throws ReadsValidationException {
        File f1 = Paths.get("src/test/resources/rawreads/EP0_GTTCCTT_S1.txt.gz").toFile();
        File f2 = Paths.get("src/test/resources/rawreads/EP0_GTTCCTT_S2.txt.gz").toFile();

        ValidatorWrapper vw = new ValidatorWrapper(Arrays.asList(f1, f2), FileFormat.FASTQ, READ_COUNT_LIMIT);
        vw.run();
    }

    @Test
    public void
    testValidPair2() throws ReadsValidationException {
        File f1 = Paths.get("src/test/resources/rawreads/CI26.26-min.fastq.gz").toFile();
        File f2 = Paths.get("src/test/resources/rawreads/RRCI26.26-min.fastq.gz").toFile();

        ValidatorWrapper vw = new ValidatorWrapper(Arrays.asList(f1, f2), FileFormat.FASTQ, READ_COUNT_LIMIT);
        vw.run();
    }

    @Test
    public void
    testValidPair3() throws IOException, ReadsValidationException {
        File output_dir = createOutputFolder();
        Path f1 = saveRandomized(
                "@A00953:544:HMTFHDSX3:2:1101:6768:1\n" +
                        "GAATAAAGATAGAGATTTGAGTAAGTATGCTCAATATAAATATTATAACCATGATATTAATGCAGGACTAAATATTATCCGCGAGAAATATAGATTGAACTTTGGAGTATCTTTGCAACCGCAGAACACCAGATTGGATTATAAGAAGGCC\n" +
                        "+\n" +
                        "FFFFFFFF:FFFFFFFFF:FFFFFFFF,FFFFFFFF,FFFFFFFFF:FFFFFFFFFFFF:F,FFF,FFFFFFFF:FFFF:FFFFFFFFFF:F:FF:FFF:::,FFFFFFFFF,:F,FFFFFFF:FFFFFFFF:FFFFFFFFFFF,FFF:FF\n" +
                        "@A00953:544:HMTFHDSX3:2:1101:16278:2\n" +
                        "CTTTTGAAAAGGATGAAATCGATGTGGACCTGGAATAAATTCAAAAAACTGATTTTTGCGGCGGGTCTGTGTTTTGCTCTTCCGATGCCGGCGCAGGCGTGCATTCCTCTGAACCCGATTTGTCTGTTCCAGATGATTCTGAAGGTTACGC\n" +
                        "+\n" +
                        "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFF\n" +
                        "@A00953:544:HMTFHDSX3:2:1101:23981:3\n" +
                        "ATCCTACACCCACGGCAGATAGGGACCGAACTGTCTCACGACGTTCTAAACCCAGCTCGCGTACCACTTTAATCGGCGAACAGCCGAACCCTTGGGACCTTCTCCAGCCCCAGGATGTGATGAGCCGACATCGAGGTGCCAAACTCCGCCG\n" +
                        "+\n" +
                        "FFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF",
                output_dir.toPath(), true, "fastq-10-7-1", "gz");
        Path f2 = saveRandomized(
                "@A00953:544:HMTFHDSX3:2:1101:6768:1\n" +
                        "GTAAAGCGCAACTGGCTGACTTTAGAAAATCTGTAACGAAAATCAACGTTAGGAGCAAAATTGAACACATTTCTTTTTACCACAGTATCTACTTCGGCCTTCTTATAATCCAATCTGGTGTTCTGCGGTTGCAAAGATACTCCAAAGTTCA\n" +
                        "+\n" +
                        ",FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFF:FFFFFFFFFFFFFF:F:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFF\n" +
                        "@A00953:544:HMTFHDSX3:2:1101:16278:2\n" +
                        "CCGGAATATGAGGAATGACGGCCGGAATGGAAACGAAGTCCATAACGGGCATTCCCGGCGTAACCTTCAGAATCATCTGGAACAGACAAATCGGGTTCAGAGGAATGCACGCCTGCGCCGGCATCGGAAGAGCAAAACACAGACCCGCCGC\n" +
                        "+\n" +
                        "FFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFF:FFFFF:FFFFFFFFFFFFFFFF\n" +
                        "@A00953:544:HMTFHDSX3:2:1101:23981:3\n" +
                        "GCATAAGCCAGCCTGACTGCGAGAGCGACAGTTCGAGCAGAGACGAAAGTCGGTCATAGTGATCCGGTGGTCCCGAGTGGAAGGGCCATCGCTCAACGGATAAAAGGTACTCCGGGGATAACAGGCTGATTCCGCCCAAGAGTTCACATCG\n" +
                        "+\n" +
                        "FF:FFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF::FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF",
                output_dir.toPath(), true, "fastq-10-7-2", "gz");

        ValidatorWrapper vw = new ValidatorWrapper(
                Arrays.asList(f1.toFile(), f2.toFile()), FileFormat.FASTQ, READ_COUNT_LIMIT);
        vw.run();
    }

    @Test
    public void
    testValidPair4() throws IOException, ReadsValidationException {
        File output_dir = createOutputFolder();
        Path f1 = saveRandomized(
                "@A00953:544:HMTFHDSX3:2:1101:6768:1814\n" +
                        "GAATAAAGATAGAGATTTGAGTAAGTATGCTCAATATAAATATTATAACCATGATATTAATGCAGGACTAAATATTATCCGCGAGAAATATAGATTGAACTTTGGAGTATCTTTGCAACCGCAGAACACCAGATTGGATTATAAGAAGGCC\n" +
                        "+\n" +
                        "FFFFFFFF:FFFFFFFFF:FFFFFFFF,FFFFFFFF,FFFFFFFFF:FFFFFFFFFFFF:F,FFF,FFFFFFFF:FFFF:FFFFFFFFFF:F:FF:FFF:::,FFFFFFFFF,:F,FFFFFFF:FFFFFFFF:FFFFFFFFFFF,FFF:FF\n" +
                        "@A00953:544:HMTFHDSX3:2:1101:16278:1815\n" +
                        "CTTTTGAAAAGGATGAAATCGATGTGGACCTGGAATAAATTCAAAAAACTGATTTTTGCGGCGGGTCTGTGTTTTGCTCTTCCGATGCCGGCGCAGGCGTGCATTCCTCTGAACCCGATTTGTCTGTTCCAGATGATTCTGAAGGTTACGC\n" +
                        "+\n" +
                        "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFF\n" +
                        "@A00953:544:HMTFHDSX3:2:1101:23981:1816\n" +
                        "ATCCTACACCCACGGCAGATAGGGACCGAACTGTCTCACGACGTTCTAAACCCAGCTCGCGTACCACTTTAATCGGCGAACAGCCGAACCCTTGGGACCTTCTCCAGCCCCAGGATGTGATGAGCCGACATCGAGGTGCCAAACTCCGCCG\n" +
                        "+\n" +
                        "FFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF",
                output_dir.toPath(), true, "fastq-10-6-1", "gz");
        Path f2 = saveRandomized(
                "@A00953:544:HMTFHDSX3:2:1101:6768:1814\n" +
                        "GTAAAGCGCAACTGGCTGACTTTAGAAAATCTGTAACGAAAATCAACGTTAGGAGCAAAATTGAACACATTTCTTTTTACCACAGTATCTACTTCGGCCTTCTTATAATCCAATCTGGTGTTCTGCGGTTGCAAAGATACTCCAAAGTTCA\n" +
                        "+\n" +
                        ",FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFF:FFFFFFFFFFFFFF:F:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFF\n" +
                        "@A00953:544:HMTFHDSX3:2:1101:16278:1815\n" +
                        "CCGGAATATGAGGAATGACGGCCGGAATGGAAACGAAGTCCATAACGGGCATTCCCGGCGTAACCTTCAGAATCATCTGGAACAGACAAATCGGGTTCAGAGGAATGCACGCCTGCGCCGGCATCGGAAGAGCAAAACACAGACCCGCCGC\n" +
                        "+\n" +
                        "FFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFF:FFFFF:FFFFFFFFFFFFFFFF\n" +
                        "@A00953:544:HMTFHDSX3:2:1101:23981:1816\n" +
                        "GCATAAGCCAGCCTGACTGCGAGAGCGACAGTTCGAGCAGAGACGAAAGTCGGTCATAGTGATCCGGTGGTCCCGAGTGGAAGGGCCATCGCTCAACGGATAAAAGGTACTCCGGGGATAACAGGCTGATTCCGCCCAAGAGTTCACATCG\n" +
                        "+\n" +
                        "FF:FFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF::FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF",
                output_dir.toPath(), true, "fastq-10-6-2", "gz");

        ValidatorWrapper vw = new ValidatorWrapper(
                Arrays.asList(f1.toFile(), f2.toFile()), FileFormat.FASTQ, READ_COUNT_LIMIT);
        vw.run();
    }

    @Test
    public void
    testValidPair5() throws IOException, ReadsValidationException {
        File output_dir = createOutputFolder();
        Path f1 = saveRandomized(
                "@JAXVAFT.MTP3.D21.48_239647 M00990:616:000000000-JKYVV:1:1101:20204:2129 1:N:0:TTACTGTGCG+GATTCCTA\n" +
                        "TACGTAGGTGGCGAGCGTTATCCGGAATGATTGGGCGTAAAGGGTGCGCAGGCGGTCCTGCAAGTCTGGAGTGAAACGCATGAGCTCAACTCATGCATGGCTTTGGAAACTGGAGGACTGGAGAGCAGGAGAGGGCGGTGGAACTCCATGTGTAGCGGTAAAATGCGTAGATATATGGAAGAACACCAGTGGCGAAGGTGGTTGCCTGGCCTGCTGCTGACGCCGAGGCACTAAAGCTGCGGGAGCAAAA\n" +
                        "+\n" +
                        "ABBBA4CABDAFG2AEEGEGFDFGGEEGGDFFHGHGAAEEADDE1BFEGGC?1E?>E@GEHHFHDHHFHBA3?344B/E//?3?DGFGHBG3??FDGHGFHBBG?C11DGDBB00//GFGHB01/<<0C<.<.AA<<--:.;CG00;CHH0CC0;0@B;C?BFFFGGGG@999C;0FBB0/:F.FF?E.9/BE--9-.9999..99BFE?F/..;:FFBB9AD----9@E../;//;///9----;A//.\n" +
                        "@JAXVAFT.MTP3.D21.48_239648 M00990:616:000000000-JKYVV:1:1101:15017:2176 1:N:0:TTACTGTGCG+GATTCCTA\n" +
                        "TACGTAGGGGGCAAGCGTTATCCGGATTTACTGGGTGTAAAGGGAGCGTAGACGGCTGAGCAAGTCAGAAGTGAAAGCCCGGGGCTCAACCCCGGGACTGCTTTTGAAACTGCGGAGCTGGAGAGCAGGAGAGGGCGGTGGAACTCCATGTGTAGCGGTAAAATGCGTAGATATATGGAAGAACACCAGTGGCGAAGGCGGCCGCCTGGCCTGTTGCTGACGCTGAGGCACGAAAGCGTGGGGAGCAAAA\n" +
                        "+\n" +
                        "ABBBBFFCCCCCGGGGGGGGGGHGGGFGHHHHGHHHGGHGHHHGHGGGFGGGHGEGEFHHGHHHHHGHHGHHHHHHGHHHGGGGGGHHHHHGGGGGCGHHHHHHHHGHHHHHHGGCGGHHHGHHGHEGHGGHGFFGGGGGGGGGGGGGGGGGGGGGDGGGGFFFFFFFDFFFFFFFFFFFFFFFFFFFFFFFFFFFFBFFFA@--CD.AFFFFFFEFFFFF?FADDFEEFFFDFFFFF:ADFFEFF.BB/\n" +
                        "@JAXVAFT.MTP3.D21.48_239649 M00990:616:000000000-JKYVV:1:1101:11786:2220 1:N:0:TTACTGTGCG+GATTCCTA\n" +
                        "TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGGATGCCAAGTCAGCGGTAAAAAAGCGGTGCTCAACGCCGTCGAGCCGTTGAAACTGGCGTTCTTGAGTGGGCGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACTCCGATTGCGAAGGCAGCATACCGGCTCCCTACTGACGCTGAGGCACGAAAGCGCGGGGAGCAAACA\n" +
                        "+\n" +
                        "1>AA?AA?AFFBE1E0AECGAGHCGGGEHG2FF1/EEEEFFHC/F/?AFGGE1E/E/EFFE1BG>FF1>EGEAF1GF?E0/@/CCGHF1?B/@/A@//-.>->.<<DD0==GE<<@CEH0CGCCACA@-9?BBBFG0F?--@GF0;.C?E--/9F?;-9-/B9BFFBFFF:FFFFF-----;9BBFB@BB-BB-@@A-F-AFFFF/--9-;E-BF-//B/AF?AAFEE-9--9-9/9-;=;;=-A-FFF9",
                output_dir.toPath(), true, "fastq-10-3-1", "gz");
        Path f2 = saveRandomized(
                "@JAXVAFT.MTP3.D21.48_239647 M00990:616:000000000-JKYVV:1:1101:20204:2129 2:N:0:TTACTGTGCG+GATTCCTA\n" +
                        "GAAGGTGACGAAGTTTGTACGGAATGATTGGACGTAAGGGAAGCGCAGACGGTCCTGCAAGTCTGGAGTGAAACGGATGAGCTGAACTCATGCATAGCTTTGGAAACTGGAAGACTAGAGAGCAGGAGAAGGCGGTGGAACTCCATATGTAGCGGTAAAAGGCGTAGATATATGGAAGAACACCAATGACGAAGGCGGCCGCCTGGCCTGTTGCTGACGCTGAGGCACGAAAGCGTGGGGAGCAAGTAAG\n" +
                        "+\n" +
                        ";///--;--///-//--9;///FB//;/-9---//;/-/..-9:0....9..;.900009:0///0:GBGFGF?<0CGC=0000<000GGEHFD=0/<1111FFFG1F1D>11<>11G1GG?111?11///CCHGGGD0BB11DFB1B//<@<B1HHF<</E/1F22BEFFGBDFGF2/F>1B220//0FEEEE@//EA/CGFA///FGHFHE/E0E01FB11EFA0EAC0E1EAB11@@1113B1>>11\n" +
                        "@JAXVAFT.MTP3.D21.48_239648 M00990:616:000000000-JKYVV:1:1101:15017:2176 2:N:0:TTACTGTGCG+GATTCCTA\n" +
                        "GTAGGGGGCCAGGGTTATCCGGATTTACTGGGTTTAAAGGGAGCGTAGACTGCTGAGCAAGTCAGAAGTGAAAGCCCGGGGCTCAACCCCGGGACTGCTTTTGAAACTGCGGAGCTGGAGAGCAGGAGAGGGCGGTGGAACTCCATGTGTAGCGGTAAAATGCGTAGATATATGGAAGAACACCAGTGGCGAAGGCGGCCGCCTGGCCTGTTGCTGACGCTGAGGCACGAAAGCGTGGGGAGCAAATAGG\n" +
                        "+\n" +
                        "..9=9999....///D;.9./9//B;/99./9///9/.E;;.-.//;99./FFFF9;/////EFFB;B/FFEGC?-?.EGFF9-;?B?GGFC/C//GGDGDGHGBE?DGDCFHF11FFG=1BHFFFC<CC/DDHGBBBBHHHFFGBHHCGFGGHFHHGGE??FGGDF4FDHFFF4FHHEGG3GHHGGGGHHGGGGE>EAEHHEHGFECGHHHHFGGFGHGHHGHGDFGGFEECEECEFFFDFFDFAABBA\n" +
                        "@JAXVAFT.MTP3.D21.48_239649 M00990:616:000000000-JKYVV:1:1101:11786:2220 2:N:0:TTACTGTGCG+GATTCCTA\n" +
                        "GTGCTGAAGGGAGGCTTATACGAGTATAATGGGTGTAAAGGGAGCGTAGGCGGGATGCCAAGGCCGCGGTAAAAAAGCGGGGTTCACGGCCGTCGGGCCGTTGAAACTGGCGTTCTTGAGTGGGCGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATATATATCACGCAGAACTCCGATTACGAAGGCAGCATACCGGCGCCCTACTGACGCTGAGGCACGACAGCGTGGGGAGCAAACAGG\n" +
                        "+\n" +
                        "-/9//-/--9;//////;---9/////9/--BB/9;//-------9---9-;////9///-9--999/9--999--999..0..-------:-----<.=0DD<=0-<-.<0F1<1000<//C?/DF>2@1//</022//CB//<00F<E///F2>F>2222F2FB222GE?/>11B1BE?CEAGF>B0>FF@CB122@EAAECEEE0E1GAB/A000B1FGAGAEAFGA0EE?AG?E11FDBDFA>>1>",
                output_dir.toPath(), true, "fastq-10-3-2", "gz");

        ValidatorWrapper vw = new ValidatorWrapper(
                Arrays.asList(f1.toFile(), f2.toFile()), FileFormat.FASTQ, READ_COUNT_LIMIT);
        vw.run();
    }

    @Test
    public void
    testValidPair6() throws IOException, ReadsValidationException {
        File output_dir = createOutputFolder();
        Path f1 = saveRandomized(
                "@JAXVAFT/1\n" +
                        "TACGTAGGTGGCGAGCGTTATCCGGAATGATTGGGCGTAAAGGGTGCGCAGGCGGTCCTGCAAGTCTGGAGTGAAACGCATGAGCTCAACTCATGCATGGCTTTGGAAACTGGAGGACTGGAGAGCAGGAGAGGGCGGTGGAACTCCATGTGTAGCGGTAAAATGCGTAGATATATGGAAGAACACCAGTGGCGAAGGTGGTTGCCTGGCCTGCTGCTGACGCCGAGGCACTAAAGCTGCGGGAGCAAAA\n" +
                        "+\n" +
                        "ABBBA4CABDAFG2AEEGEGFDFGGEEGGDFFHGHGAAEEADDE1BFEGGC?1E?>E@GEHHFHDHHFHBA3?344B/E//?3?DGFGHBG3??FDGHGFHBBG?C11DGDBB00//GFGHB01/<<0C<.<.AA<<--:.;CG00;CHH0CC0;0@B;C?BFFFGGGG@999C;0FBB0/:F.FF?E.9/BE--9-.9999..99BFE?F/..;:FFBB9AD----9@E../;//;///9----;A//.",
                output_dir.toPath(), true, "fastq-10-4-1", "gz");
        Path f2 = saveRandomized(
                "@JAXVAFT/2\n" +
                        "GAAGGTGACGAAGTTTGTACGGAATGATTGGACGTAAGGGAAGCGCAGACGGTCCTGCAAGTCTGGAGTGAAACGGATGAGCTGAACTCATGCATAGCTTTGGAAACTGGAAGACTAGAGAGCAGGAGAAGGCGGTGGAACTCCATATGTAGCGGTAAAAGGCGTAGATATATGGAAGAACACCAATGACGAAGGCGGCCGCCTGGCCTGTTGCTGACGCTGAGGCACGAAAGCGTGGGGAGCAAGTAAG\n" +
                        "+\n" +
                        ";///--;--///-//--9;///FB//;/-9---//;/-/..-9:0....9..;.900009:0///0:GBGFGF?<0CGC=0000<000GGEHFD=0/<1111FFFG1F1D>11<>11G1GG?111?11///CCHGGGD0BB11DFB1B//<@<B1HHF<</E/1F22BEFFGBDFGF2/F>1B220//0FEEEE@//EA/CGFA///FGHFHE/E0E01FB11EFA0EAC0E1EAB11@@1113B1>>11",
                output_dir.toPath(), true, "fastq-10-4-2", "gz");

        ValidatorWrapper vw = new ValidatorWrapper(
                Arrays.asList(f1.toFile(), f2.toFile()), FileFormat.FASTQ, READ_COUNT_LIMIT);
        vw.run();
    }

    @Test
    public void
    testValidPair7() throws IOException, ReadsValidationException {
        File output_dir = createOutputFolder();
        Path f1 = saveRandomized(
                "@E00528:414:HVNJLCCXY:1:1101:7598:1854 1:N:0\n" +
                        "CTTTCAGGTTGTTGAGAATGTCGTTAGCTTCTTGTAAGGTATTGCGTCCCTTTTTAGCAGCTTCTTCAGCAAGAGCTTTGGCAGCGTCGGCTCGAGCTAGGAGTTGGTCAGCGGTCTGCTGTTCGGCTTTCCCCTTCTCTAGAAGGTTCT\n" +
                        "+\n" +
                        "```e``Vi`ei`eiLiVLV[V[L[Liiiiiiiii`eVeieeeieeiieiieiei`iiiiiiiiieiieiieiieii[eieiiiiiiiiiiiiiiiiie[LeeeL`iiiie[iLeeiiiiiieeLeiii[ieLeieiiiLee[iVeii[L[\n" +
                        "@E00528:414:HVNJLCCXY:1:1101:12773:1854 1:N:0\n" +
                        "CAAGGAACCTCTGTGGCGTTCACAAAAAACCCCGCACCAGTGCCGAAGTCCCAGCTGTCATCTTCTCCCTTATTATAGCAGCCACGGGGGCTGGTATCAGGAGCAATGACCATAAGGCCATGTTCTGAGGCAGCTTGTTGACAGCCAGAC\n" +
                        "+\n" +
                        "``L``eLiiiiVeeeiiee[[eVeVLeVeieiiiL`ViiieiiiiieieeVeL[[`[ViieiiiiiiiL[eeL[[eLV`ieieL`ii[ieL`iiiee`iiL`VL`e``[V[LLV`Lie`[e`i[`e`ieiiL[[iiii`ie``L[`Le`e\n" +
                        "@E00528:414:HVNJLCCXY:1:1101:17421:1889 1:N:0\n" +
                        "TTCATGGCTAGGACATTTTCTAAACAGTAAAACGACTTCAGGATATATCTCTTCGACATTAAATGTAGCAAAATGGAAAAAAAGATTTACAAGTAGATAAAGTAATAACAGAATATTACAAAGCAGAGAACACACACACAAAAGAATCAA\n" +
                        "+\n" +
                        "``LeeiL`Ve``ieeii[`Leeei[iiiiiiiiiiiLeieiiiiiiie`iiiiiiiiiiiiiiiiiiiiiii`Li``eiiiiiii`e[eeeiiiiii`LeViL`iiiieiiiiiiiiiiiiiiiiiee[eieeieiiiVeiiieee`iee",
                output_dir.toPath(), true, "fastq-10-1-1", "gz");
        Path f2 = saveRandomized(
                "@E00528:414:HVNJLCCXY:1:1101:7598:1854 2:N:0\n" +
                        "ATAAAATCAAGAAAGAAGCTGCAGATCTGGACCGTCTGATTGAACAGAAGCTAAAAGATTAGGAGGACCTTAGGGAAGACATGCGAGGAAAGGAACATGAAGTAAAGAACGTTCTAGAGAAGGGGAAAGCCGAATAGCAGATCGCTGACC\n" +
                        "+\n" +
                        "``eeei[iieeiieieiiV[[eiii`L[eeieiieiVeiLeiiLLee[eiVeeiiiiiiiLL[eiiiiVVV`[V``iiiVeV`LVeiiiiiiiVeiiii`eeVV[`ie[[LVVVVe`[`eiiiV`V`[[L``[iLV[`[e`LV`V`V`[e\n" +
                        "@E00528:414:HVNJLCCXY:1:1101:12773:1854 2:N:0\n" +
                        "CCTGCACTTTACTGGCTGTCTGGTTTAAATTGCACAGAACAAAATTTCATATCAAAGTCTGGCTGTCAACAAGATGCCTTAGAACATGGCCTTGTGGTGATTGCTCCTGATACCAGCAGCCGTGGCTGCAATATTAAAGGAGAAGATGAC\n" +
                        "+\n" +
                        "```e`iL`eiie`iiL`eVL[eiiii`eL[ee[eiieii[eie`ei`iiiLi`ieeeeiii[ei``eiLVL`eLL[`iiLV`Vieiii`iiii[eiieLL``eL`eV`ieV`eL[eLVLV`i[``HL[HVV[eeieiLV[[VeeVV`eie\n" +
                        "@E00528:414:HVNJLCCXY:1:1101:17421:1889 2:N:0\n" +
                        "CAGAATGCTTAAAGGCCGTCTGTAGGAGGCTCCAGTACATGGAAAGAAAAGGATTCAACACAGTGTGGTCATATGATAAATAAGTGATTTATAAACAAACAAGAGTGATATTTTGCTAGTTAACAAAGTTAACAGATCTTCTAGCCTCAT\n" +
                        "+\n" +
                        "```eeiLeieiiiii`eiiL[eeiiiL[LVe`Veie`iiLeieLVeeiiiiii`ieiiieeL[[e`iiV`VL`Veeiiieeeiiiiiie`iiiLLiiieeLeeii`eee`Veieee[L``iiiVeL``e`ie`VVeLLV`[V``[eLHL[",
                output_dir.toPath(), true, "fastq-10-1-2", "gz");

        ValidatorWrapper vw = new ValidatorWrapper(
                Arrays.asList(f1.toFile(), f2.toFile()), FileFormat.FASTQ, READ_COUNT_LIMIT);
        vw.run();
    }

    @Test
    public void
    testInvalidPair() throws IOException {
        File output_dir = createOutputFolder();
        Path f1 = saveRandomized("@NAME1\nACGT\n+\n1234\n"
                + "@NAME1\nACGT\n+\n1234", output_dir.toPath(), true, "fastq-1", "gz");
        Path f2 = saveRandomized(
                "@NAME2\nACGT\n+\n1234", output_dir.toPath(), true, "fastq-2", "gz");

        try {
            ValidatorWrapper vw = new ValidatorWrapper(
                    Arrays.asList(f1.toFile(), f2.toFile()), FileFormat.FASTQ, READ_COUNT_LIMIT);
            vw.run();
            fail();
        } catch (ReadsValidationException e) {
            e.printStackTrace();
        }
    }

    @Test
    public void
    testInvalidPair2() throws IOException {
        File output_dir = createOutputFolder();
        Path f1 = saveRandomized("@NAME1\nACGT\n+\n1234", output_dir.toPath(), true, "fastq-1", "gz");
        Path f2 = saveRandomized("@NAME2\nACGT\n+\n1234\n"
                + "@NAME2\nACGT\n+\n2341", output_dir.toPath(), true, "fastq-2", "gz");

        try {
            ValidatorWrapper vw = new ValidatorWrapper(
                    Arrays.asList(f1.toFile(), f2.toFile()), FileFormat.FASTQ, READ_COUNT_LIMIT);
            vw.run();
            fail();
        } catch (ReadsValidationException e) {
            e.printStackTrace();
        }
    }

    @Test
    public void
    testInvalidPair3() throws IOException {
        File output_dir = createOutputFolder();
        Path f1 = saveRandomized("@NAME\nACGT\n+\n1234", output_dir.toPath(), true, "fastq-1", "gz");
        Path f2 = saveRandomized("@NAME/2\nACGT\n+\n1234\n"
                + "@NAME/3\nACGT\n+\n2341", output_dir.toPath(), true, "fastq-2", "gz");

        try {
            ValidatorWrapper vw = new ValidatorWrapper(
                    Arrays.asList(f1.toFile(), f2.toFile()), FileFormat.FASTQ, READ_COUNT_LIMIT);
            vw.run();
            fail();
        } catch (ReadsValidationException e) {
            e.printStackTrace();
        }
    }

    @Test
    public void
    testInvalidPair4() throws IOException {
        File output_dir = createOutputFolder();
        Path f1 = saveRandomized("@NAME1/1\nACGT\n+\n1234", output_dir.toPath(), true, "fastq-1", "gz");
        Path f2 = saveRandomized("@NAME2/2\nACGT\n+\n1234\n"
                + "@NAME/2\nACGT\n+\n2341", output_dir.toPath(), true, "fastq-2", "gz");

        try {
            ValidatorWrapper vw = new ValidatorWrapper(
                    Arrays.asList(f1.toFile(), f2.toFile()), FileFormat.FASTQ, READ_COUNT_LIMIT);
            vw.run();
            fail();
        } catch (ReadsValidationException e) {
            e.printStackTrace();
        }
    }

    @Test
    public void
    testInvalidPair5() throws IOException {
        File output_dir = createOutputFolder();
        Path f1 = saveRandomized(
                "@FCC4696:1:1101:19283:1592#TTC_CGTGACGG_CTAGTTAT/1\n" +
                        "ACTCCTACGGGAGGCAGCAGTGAGGAATATTGGTCAATGGGCGCTAGCCTGAACCAGCCAAGTAGCGTGAAGGATGACTGCCCTATGGGTTGTAAACTTCTTTTATATGGGAATAAAGTGCAGTATGTATACTGCTTTGCATGTACCTTATGAATAAGGATCGGCTAACTCCGTTCCAGCAGCCGCGGTCATACGTAGGATCCGAGCGTTATCCGGATTTATTTGGTTTTAATGGATCGT\n" +
                        "+\n" +
                        "CCGGGGGGGGGGGGGGGGGGGGGGGGGFFGGF9FCA,EFF,6+@F76C<F8,<ED,,:C,,,C,,CCGGC7@D@FE8@C,CEEFDE,,<FCCF,,,BFEFEGGGAF9@,5=,,C,,,,C,BE,A,C<F9F,EE<?FFE5?,A9F,4AFF9F,,,:,,,9,6@++@D+,@ADFCF8DB9,@,,6=8C=8@+4@,@8,,5*,7<1***31C9*;:B***2=D+=C8++?CG*2*+3+*+/28\n" +
                        "@FCC4696:1:1101:8646:1667#TTC_CGTGACGG_CTAGTTAT/1\n" +
                        "ACTCCTACGGGAGGCAGCAGTGAGGAATATTGGTCAATGGGCGGGAGCCTGAACCAGCCAAGTCGCGTGAGGGATGACGGCCCTATGGGTTGTAAACCTCTTTTGCCGGGGAGCAAAGTGCCGCACGTGTGCGGTTTGGAGAGTACCCGGAGAAAAAGCATCGGCTAACTCCGTTCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTTGGTTTAAAGGGTGCGTAGGCGGACGCTTAAGTCAGCGGTAAAATTTCGGGGCTCAACCTCGTCTAGCCGTTGA\n" +
                        "+\n" +
                        "CCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGF<7:B7CFFGE@8FG8,CF8CFGF@GGGGGGFGEGGFF7:CFFGGGCGGGGGGA9<FFGGFGGGGGGG7CG7=FD8,,FBFGD@+>@FFF+@CFGGG>7,,>CF9CCFECFC*,3,<<C,?<;:BF*+?>>FCGFGF9CE**<;:EEGC*;F7E=***85?6C)*99:DC3CDGDC*9CD7FD*2:CDG4*0C457)-9?69:7(((-2;3:(4<D464:29:(14)6:--2(,3(414>:3311<?(,6>?0:<,2\n" +
                        "@FCC4696:1:1101:9652:1672#TTC_CGTGACGG_TTAGTTAT/1\n" +
                        "ACTCCTACGGGAGGCAGCAGTGAGGAATATTGGTCAATGGGCGAGAGCCTGAACCAGCCAAGTCGCGTGAGGGAGTACTGCCCTATGGGTTGTAAACCTCTTTTGTCGGGGAGCAAAAGCCGGACGTGTCCGGCTGTGAGAGTACCCGAAGAAAAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTCATACGGAGGATGCGAGCGTTATCCGGATTTATTGTGTTTAAAGGGTGCGCAGGCGGATTTTTAAGTCAGCGGTCAAATACG\n" +
                        "+\n" +
                        "CCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGF7CCCEGEGC@FGGC,CDCCFGFDGGGGGGGGGGGGGFCEFFGGGGGGGGGGDFGGGGGGGGGGGGGDGGCGGG,,EBFGG7:CEGGGGFECFGEG,,<=FFFGEG>7:F@7,@@CGGG@*3BE88DFGGGGF:B9>FFCBFGG=GG*3C9CEC5:=EGGG553DDGGDGFGD33CDF7DF7*:>GF>*0*:B>577;(43;(1(4:?<>:7>BA))519?1()6A:6,\n",
                output_dir.toPath(), true, "fastq-1", "gz");
        Path f2 = saveRandomized(
                "@FCC4696:1:1101:19283:1592#C_CGTGACGG_CTAGTTAT/2\n" +
                        "GGACTACAGGGTTATCTAATCCTGTTTGATACCCACACTTTCGAGCATCAGCGTCAGTTACAGTCCAGCAGGCTGCCTTCGCAATCGGATTTCTTCTTGATATCTAAGCATTTCACCGCTACACCACGAATTCCGCCTACCTCTCCTGTACTCAACGCTGCCAGTATCTACTGCAATTTTACGGTTGCGCCGCACACTTTCCCACCTGACTTAACAACCCGCCTACGCTCCCTTTACACC\n" +
                        "+\n" +
                        ",,,AF8F9,,,,;,CCC,,CE@E,CEF8,C,CEE9F,BFFGD7,,8@FFE9<+CB7,FF,C9@FFE,,6,96:?,?AFGG,?,+@E@,,,8FBFGG,@,,C,9=@,,,@,@EED,@E76@+@,@B+@7++9=EF?=8A+ADFFD*;D?@8DDD+**3*60=D=+@7?;++4<;<++3:;C*;*6ADC++00):))10:AC@*00)/:))7A<**1*157225>?2:227?91@93)))/:\n" +
                        "@FCC4696:1:1101:8646:1667#C_CGTGACGG_CTAGTTAT/2\n" +
                        "GGACTACCCGGGTATCTAATCCTGTTCGATACCCGCACTTTCGTGCCTCAGCGTCAGTTGAGCGCCGGTATGCTGCCTTCGCAATCGGGGTTCTGCGTGATATCTATGCATTTCACCGCTACACCACGCATTCCGCGTACTTCTCGCCCACTCAAGGCACCCAGTTTCAACGGCTCGACGAGGTTGAGCCCCGCAATTTTACCGCTGACTTACGCGTCCGCCTACGCACCCTTTCACCCCAATAAATCCGGATAACGCTCTCATCCTCCG\n" +
                        "+\n" +
                        "A<ACGGGGGCCFGGGGGCCGGGGAFGGGGFFGGGGGDGGGGGGG<FGGGGGCEGGGGGG9ECGGGGFCFGGGGG,EFGGGGGC7FGFD+:CGGGEDCE@?FFGEEFFEFFFGGGFGGCBCBCDGF?DEE>CGGGDEGECGGGFFF>AB?CEAF?5**39?7?CGFFG<+1=:49:DF(-5(:8CFGGGCF53)3(,6A?EC<3603967AA<)(324>>>BF?,:<9(79<:B7)6).8>?<::)4:62,(-2(98:>(-(,-6<A?-42\n" +
                        "@FCC4696:1:1101:9652:1672#C_CGTGACGG_TTAGTTAT/2\n" +
                        "GGACTACTCGGGTATCTAATCCTGTTCGATACCCACGCTTTCGTGCTTCAGCGTCAGTTGGGCGCCGGTATGCTGCCTTCGCAATCGGAGTTCTGCGTGATATCTATGCATTTCACCGCTACACCACGCATTCCGCGTACTTCTCGCCCACTCAAGAACGCCAGTTTCAACGGCGGACGCAGGTTGAGCCCCCGTATTTGACCGCTGACTTAAAAATCCGCCTGCGCACCCTTTTAACCC\n" +
                        "+\n" +
                        "@@CCGGGGDB@EGGGGG9EGGGGCGGGGGGGGGGGGGGGDGGGF@FGGGGFCFGGGGGGFE7FGGG>CFGGFGFFFGGGGGCFFGGFD,+FFGGC;CE7CGGGCDFGFGGEGGGGGFDBECEFFG6DEG6EGGGGGDG>EGGF9DBA?BCDFF5A**3A7EGDGGFGF6+3:??@BAB(:1DEGGGGFFGD?E><B;E4)-5,46>6>72:)2)6)3:?BB?B?6224>BB:?B))).56\n",
                output_dir.toPath(), true, "fastq-2", "gz");

        try {
            ValidatorWrapper vw = new ValidatorWrapper(
                    Arrays.asList(f1.toFile(), f2.toFile()), FileFormat.FASTQ, READ_COUNT_LIMIT);
            vw.run();
            fail();
        } catch (ReadsValidationException e) {
            e.printStackTrace();
        }
    }

    @Test
    public void
    testEmptyReadDoubleFile() throws IOException {
        File output_dir = createOutputFolder();
        Path f1 = saveRandomized(
                "@E00528:414:HVNJLCCXY:1:1101:7598:1854 1:N:0\n" +
                        "CTTTCAGGTTGTTGAGAATGTCGTTAGCTTCTTGTAAGGTATTGCGTCCCTTTTTAGCAGCTTCTTCAGCAAGAGCTTTGGCAGCGTCGGCTCGAGCTAGGAGTTGGTCAGCGGTCTGCTGTTCGGCTTTCCCCTTCTCTAGAAGGTTCT\n" +
                        "+\n" +
                        "```e``Vi`ei`eiLiVLV[V[L[Liiiiiiiii`eVeieeeieeiieiieiei`iiiiiiiiieiieiieiieii[eieiiiiiiiiiiiiiiiiie[LeeeL`iiiie[iLeeiiiiiieeLeiii[ieLeieiiiLee[iVeii[L[\n" +
                        "@E00528:414:HVNJLCCXY:1:1101:12773:1854 1:N:0\n" +
                        "CAAGGAACCTCTGTGGCGTTCACAAAAAACCCCGCACCAGTGCCGAAGTCCCAGCTGTCATCTTCTCCCTTATTATAGCAGCCACGGGGGCTGGTATCAGGAGCAATGACCATAAGGCCATGTTCTGAGGCAGCTTGTTGACAGCCAGAC\n" +
                        "+\n" +
                        "``L``eLiiiiVeeeiiee[[eVeVLeVeieiiiL`ViiieiiiiieieeVeL[[`[ViieiiiiiiiL[eeL[[eLV`ieieL`ii[ieL`iiiee`iiL`VL`e``[V[LLV`Lie`[e`i[`e`ieiiL[[iiii`ie``L[`Le`e\n" +
                        "@E00528:414:HVNJLCCXY:1:1101:17421:1889 1:N:0\n" +
                        "TTCATGGCTAGGACATTTTCTAAACAGTAAAACGACTTCAGGATATATCTCTTCGACATTAAATGTAGCAAAATGGAAAAAAAGATTTACAAGTAGATAAAGTAATAACAGAATATTACAAAGCAGAGAACACACACACAAAAGAATCAA\n" +
                        "+\n" +
                        "``LeeiL`Ve``ieeii[`Leeei[iiiiiiiiiiiLeieiiiiiiie`iiiiiiiiiiiiiiiiiiiiiii`Li``eiiiiiii`e[eeeiiiiii`LeViL`iiiieiiiiiiiiiiiiiiiiiee[eieeieiiiVeiiieee`iee",
                output_dir.toPath(), true, "fastq-10-1-1", "gz");
        Path f2 = saveRandomized(
                "@E00528:414:HVNJLCCXY:1:1101:7598:1854 2:N:0\n" +
                        "ATAAAATCAAGAAAGAAGCTGCAGATCTGGACCGTCTGATTGAACAGAAGCTAAAAGATTAGGAGGACCTTAGGGAAGACATGCGAGGAAAGGAACATGAAGTAAAGAACGTTCTAGAGAAGGGGAAAGCCGAATAGCAGATCGCTGACC\n" +
                        "+\n" +
                        "``eeei[iieeiieieiiV[[eiii`L[eeieiieiVeiLeiiLLee[eiVeeiiiiiiiLL[eiiiiVVV`[V``iiiVeV`LVeiiiiiiiVeiiii`eeVV[`ie[[LVVVVe`[`eiiiV`V`[[L``[iLV[`[e`LV`V`V`[e\n" +
                        "@E00528:414:HVNJLCCXY:1:1101:12773:1854 2:N:0\n" +
                        "\n" +
                        "+\n" +
                        "\n" +
                        "@E00528:414:HVNJLCCXY:1:1101:17421:1889 2:N:0\n" +
                        "CAGAATGCTTAAAGGCCGTCTGTAGGAGGCTCCAGTACATGGAAAGAAAAGGATTCAACACAGTGTGGTCATATGATAAATAAGTGATTTATAAACAAACAAGAGTGATATTTTGCTAGTTAACAAAGTTAACAGATCTTCTAGCCTCAT\n" +
                        "+\n" +
                        "```eeiLeieiiiii`eiiL[eeiiiL[LVe`Veie`iiLeieLVeeiiiiii`ieiiieeL[[e`iiV`VL`Veeiiieeeiiiiiie`iiiLLiiieeLeeii`eee`Veieee[L``iiiVeL``e`ie`VVeLLV`[V``[eLHL[",
                output_dir.toPath(), true, "fastq-10-1-2", "gz");

        try {
            ValidatorWrapper vw = new ValidatorWrapper(
                    Arrays.asList(f1.toFile(), f2.toFile()), FileFormat.FASTQ, READ_COUNT_LIMIT);
            vw.run();
            fail();
        } catch (ReadsValidationException e) {
            e.printStackTrace();
        }
    }

    @Test
    public void
    testEmptyFileDoubleFile() throws IOException {
        File output_dir = createOutputFolder();
        Path f1 = saveRandomized(
                "@E00528:414:HVNJLCCXY:1:1101:7598:1854 1:N:0\n" +
                        "CTTTCAGGTTGTTGAGAATGTCGTTAGCTTCTTGTAAGGTATTGCGTCCCTTTTTAGCAGCTTCTTCAGCAAGAGCTTTGGCAGCGTCGGCTCGAGCTAGGAGTTGGTCAGCGGTCTGCTGTTCGGCTTTCCCCTTCTCTAGAAGGTTCT\n" +
                        "+\n" +
                        "```e``Vi`ei`eiLiVLV[V[L[Liiiiiiiii`eVeieeeieeiieiieiei`iiiiiiiiieiieiieiieii[eieiiiiiiiiiiiiiiiiie[LeeeL`iiiie[iLeeiiiiiieeLeiii[ieLeieiiiLee[iVeii[L[\n" +
                        "@E00528:414:HVNJLCCXY:1:1101:12773:1854 1:N:0\n" +
                        "CAAGGAACCTCTGTGGCGTTCACAAAAAACCCCGCACCAGTGCCGAAGTCCCAGCTGTCATCTTCTCCCTTATTATAGCAGCCACGGGGGCTGGTATCAGGAGCAATGACCATAAGGCCATGTTCTGAGGCAGCTTGTTGACAGCCAGAC\n" +
                        "+\n" +
                        "``L``eLiiiiVeeeiiee[[eVeVLeVeieiiiL`ViiieiiiiieieeVeL[[`[ViieiiiiiiiL[eeL[[eLV`ieieL`ii[ieL`iiiee`iiL`VL`e``[V[LLV`Lie`[e`i[`e`ieiiL[[iiii`ie``L[`Le`e\n" +
                        "@E00528:414:HVNJLCCXY:1:1101:17421:1889 1:N:0\n" +
                        "TTCATGGCTAGGACATTTTCTAAACAGTAAAACGACTTCAGGATATATCTCTTCGACATTAAATGTAGCAAAATGGAAAAAAAGATTTACAAGTAGATAAAGTAATAACAGAATATTACAAAGCAGAGAACACACACACAAAAGAATCAA\n" +
                        "+\n" +
                        "``LeeiL`Ve``ieeii[`Leeei[iiiiiiiiiiiLeieiiiiiiie`iiiiiiiiiiiiiiiiiiiiiii`Li``eiiiiiii`e[eeeiiiiii`LeViL`iiiieiiiiiiiiiiiiiiiiiee[eieeieiiiVeiiieee`iee",
                output_dir.toPath(), true, "fastq-10-1-1", "gz");
        Path f2 = saveRandomized(
                "",
                output_dir.toPath(), true, "fastq-10-1-2", "gz");

        try {
            ValidatorWrapper vw = new ValidatorWrapper(
                    Arrays.asList(f1.toFile(), f2.toFile()), FileFormat.FASTQ, READ_COUNT_LIMIT);
            vw.run();
            fail();
        } catch (ReadsValidationException e) {
            e.printStackTrace();
        }
    }

    @Test
    public void
    testMultiFileValid() throws ReadsValidationException {
        Path f1 = Paths.get("src/test/resources/10x/4fastq/I1.fastq");
        Path f2 = Paths.get("src/test/resources/10x/4fastq/R1.fastq");
        Path f3 = Paths.get("src/test/resources/10x/4fastq/R2.fastq");
        Path f4 = Paths.get("src/test/resources/10x/4fastq/R3.fastq");

        ValidatorWrapper vw = new ValidatorWrapper(
                Arrays.asList(f1.toFile(), f2.toFile(), f3.toFile(), f4.toFile()), FileFormat.FASTQ, READ_COUNT_LIMIT);
        vw.run();
    }

    @Test
    public void
    testMultiFileValid2() throws IOException, ReadsValidationException {
        File output_dir = createOutputFolder();
        Path f1 = saveRandomized(
                "@A00730:546:HWCTCDRXY:2:2101:1090:1031 1:N:0:ACAGCAAC\n" +
                        "ANAGTGGGCGTGGCACGGGCGTGGCCTGCGCGCCGTGGGGCTGCGGCGGC\n" +
                        "+\n" +
                        "F#FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF",
                output_dir.toPath(), true, "fastq-10-2-1", "gz");
        Path f2 = saveRandomized(
                "@A00730:546:HWCTCDRXY:2:2101:1090:1031 2:N:0:ACAGCAAC\n" +
                        "CAGACGCGCATTAGTCTAAGGACC\n" +
                        "+\n" +
                        "FFFFFFFFFF:FFF,FFF:FFFFF",
                output_dir.toPath(), true, "fastq-10-2-2", "gz");
        Path f3 = saveRandomized(
                "@A00730:546:HWCTCDRXY:2:2101:1090:1031 3:N:0:ACAGCAAC\n" +
                        "CGCCCGCTCCTTGCGCACACGCCAGGCAGCGCCGCCGCAGCCCCACGGCG\n" +
                        "+\n" +
                        "FFFFFFFFFFFFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFFFFFF",
                output_dir.toPath(), true, "fastq-10-2-3", "gz");
        Path f4 = saveRandomized(
                "@A00730:546:HWCTCDRXY:2:2101:1090:1031 1:N:0:ACAGCAAC\n" +
                        "ACAGCAAC\n" +
                        "+\n" +
                        "FFFFFFFF",
                output_dir.toPath(), true, "fastq-10-2-4", "gz");

        ValidatorWrapper vw = new ValidatorWrapper(
                Arrays.asList(f1.toFile(), f2.toFile(), f3.toFile(), f4.toFile()),
                FileFormat.FASTQ, READ_COUNT_LIMIT);
        vw.run();
    }

    @Test
    public void
    testMultiFileValid3() throws IOException, ReadsValidationException {
        File output_dir = createOutputFolder();
        Path f1 = saveRandomized(
                "@1:1101:1027:4805 1:N:0:TTACGGTGTCATGAGG\n" +
                        "TACGTAGGTGGCGAGCGTTATCCGGAATGATTGGGCGTAAAGGGTGCGCAGGCGGTCCTGCAAGTCTGGAGTGAAACGCATGAGCTCAACTCATGCATGGCTTTGGAAACTGGAGGACTGGAGAGCAGGAGAGGGCGGTGGAACTCCATGTGTAGCGGTAAAATGCGTAGATATATGGAAGAACACCAGTGGCGAAGGTGGTTGCCTGGCCTGCTGCTGACGCCGAGGCACTAAAGCTGCGGGAGCAAAA\n" +
                        "+\n" +
                        "ABBBA4CABDAFG2AEEGEGFDFGGEEGGDFFHGHGAAEEADDE1BFEGGC?1E?>E@GEHHFHDHHFHBA3?344B/E//?3?DGFGHBG3??FDGHGFHBBG?C11DGDBB00//GFGHB01/<<0C<.<.AA<<--:.;CG00;CHH0CC0;0@B;C?BFFFGGGG@999C;0FBB0/:F.FF?E.9/BE--9-.9999..99BFE?F/..;:FFBB9AD----9@E../;//;///9----;A//.",
                output_dir.toPath(), true, "fastq-10-5-1", "gz");
        Path f2 = saveRandomized(
                "@1:1101:1027:4805 2:N:0:TTACGGTGTCATGAGG\n" +
                        "GAAGGTGACGAAGTTTGTACGGAATGATTGGACGTAAGGGAAGCGCAGACGGTCCTGCAAGTCTGGAGTGAAACGGATGAGCTGAACTCATGCATAGCTTTGGAAACTGGAAGACTAGAGAGCAGGAGAAGGCGGTGGAACTCCATATGTAGCGGTAAAAGGCGTAGATATATGGAAGAACACCAATGACGAAGGCGGCCGCCTGGCCTGTTGCTGACGCTGAGGCACGAAAGCGTGGGGAGCAAGTAAG\n" +
                        "+\n" +
                        ";///--;--///-//--9;///FB//;/-9---//;/-/..-9:0....9..;.900009:0///0:GBGFGF?<0CGC=0000<000GGEHFD=0/<1111FFFG1F1D>11<>11G1GG?111?11///CCHGGGD0BB11DFB1B//<@<B1HHF<</E/1F22BEFFGBDFGF2/F>1B220//0FEEEE@//EA/CGFA///FGHFHE/E0E01FB11EFA0EAC0E1EAB11@@1113B1>>11",
                output_dir.toPath(), true, "fastq-10-5-2", "gz");
        Path f3 = saveRandomized(
                "@1:1101:1027:4805 3:N:0:TTACGGTGTCATGAGG\n" +
                        "GAAGGTGACGAAGTTTGTACGGAATGATTGGACGTAAGGGAAGCGCAGACGGTCCTGCAAGTCTGGAGTGAAACGGATGAGCTGAACTCATGCATAGCTTTGGAAACTGGAAGACTAGAGAGCAGGAGAAGGCGGTGGAACTCCATATGTAGCGGTAAAAGGCGTAGATATATGGAAGAACACCAATGACGAAGGCGGCCGCCTGGCCTGTTGCTGACGCTGAGGCACGAAAGCGTGGGGAGCAAGTAAG\n" +
                        "+\n" +
                        ";///--;--///-//--9;///FB//;/-9---//;/-/..-9:0....9..;.900009:0///0:GBGFGF?<0CGC=0000<000GGEHFD=0/<1111FFFG1F1D>11<>11G1GG?111?11///CCHGGGD0BB11DFB1B//<@<B1HHF<</E/1F22BEFFGBDFGF2/F>1B220//0FEEEE@//EA/CGFA///FGHFHE/E0E01FB11EFA0EAC0E1EAB11@@1113B1>>11",
                output_dir.toPath(), true, "fastq-10-5-3", "gz");

        ValidatorWrapper vw = new ValidatorWrapper(
                Arrays.asList(f1.toFile(), f2.toFile(), f3.toFile()), FileFormat.FASTQ, READ_COUNT_LIMIT);
        vw.run();
    }

    @Test
    public void testMultiplePairedFastqsWithLowPairingPercentage() throws IOException {
        File output_dir = createOutputFolder();

        //Following files' pairing arrangment is:
        //f1 & f2 = 100%
        //f1 & f3 = 0%
        //f1 & f4 = 50%

        Path f1 = saveRandomized("@NAME1/1\nACGT\n+\n1234\n"
                + "@NAME2/1\nACGT\n+\n1234", output_dir.toPath(), true, "fastq-1", "gz");

        Path f2 = saveRandomized("@NAME1/2\nACGT\n+\n1234\n"
                + "@NAME2/2\nACGT\n+\n2341", output_dir.toPath(), true, "fastq-2", "gz");

        Path f3 = saveRandomized("@NAME3/3\nACGT\n+\n1234\n"
                + "@NAME4/3\nACGT\n+\n2341", output_dir.toPath(), true, "fastq-3", "gz");

        Path f4 = saveRandomized("@NAME2/4\nACGT\n+\n1234\n"
                + "@NAME4/4\nACGT\n+\n2341", output_dir.toPath(), true, "fastq-4", "gz");

        try {
            ValidatorWrapper vw = new ValidatorWrapper(
                    Arrays.asList(f1.toFile(), f2.toFile(), f3.toFile(), f4.toFile()), FileFormat.FASTQ, READ_COUNT_LIMIT);
            vw.run();
            fail();
        } catch (ReadsValidationException e) {
            e.printStackTrace();
        }
    }

    @Test
    public void testPairingThresholdPass() throws IOException, ReadsValidationException {
        File output_dir = createOutputFolder();
        Path f1 = saveRandomized(
                "@NAME1/1\nACGT\n+\n1234\n" +
                        "@NAME2/1\nACGT\n+\n1234\n" +
                        "@NAME3/1\nACGT\n+\n1234\n" +
                        "@NAME4/1\nACGT\n+\n1234\n" +
                        "@NAME5/1\nACGT\n+\n1234\n" +
                        "@NAME6/1\nACGT\n+\n1234\n",
                output_dir.toPath(), true, "fastq-1", "gz");
        Path f2 = saveRandomized(
                "@NAME1/2\nACGT\n+\n1234\n" +
                        "@NAME2/2\nACGT\n+\n1234\n" +
                        "@NAME8/2\nACGT\n+\n1234\n" +
                        "@NAME9/2\nACGT\n+\n1234\n" +
                        "@NAME10/2\nACGT\n+\n1234\n" +
                        "@NAME11/2\nACGT\n+\n1234\n",
                output_dir.toPath(), true, "fastq-2", "gz");

        ValidatorWrapper vw = new ValidatorWrapper(
                Arrays.asList(f1.toFile(), f2.toFile()), FileFormat.FASTQ, READ_COUNT_LIMIT);
        vw.run();
    }

    @Test
    public void testPairingThresholdFail() throws IOException {
        File output_dir = createOutputFolder();
        Path f1 = saveRandomized(
                "@NAME1/1\nACGT\n+\n1234\n" +
                        "@NAME2/1\nACGT\n+\n1234\n" +
                        "@NAME3/1\nACGT\n+\n1234\n" +
                        "@NAME4/1\nACGT\n+\n1234\n" +
                        "@NAME5/1\nACGT\n+\n1234\n" +
                        "@NAME6/1\nACGT\n+\n1234\n",
                output_dir.toPath(), true, "fastq-1", "gz");
        Path f2 = saveRandomized(
                "@NAME1/2\nACGT\n+\n1234\n" +
                        "@NAME7/2\nACGT\n+\n1234\n" +
                        "@NAME8/2\nACGT\n+\n1234\n" +
                        "@NAME9/2\nACGT\n+\n1234\n" +
                        "@NAME10/2\nACGT\n+\n1234\n" +
                        "@NAME11/2\nACGT\n+\n1234\n",
                output_dir.toPath(), true, "fastq-2", "gz");

        try {
            ValidatorWrapper vw = new ValidatorWrapper(
                    Arrays.asList(f1.toFile(), f2.toFile()), FileFormat.FASTQ, READ_COUNT_LIMIT);
            vw.run();
            fail();
        } catch (ReadsValidationException e) {
            e.printStackTrace();
        }
    }

    @Ignore("Manual test")
    @Test
    public void testERR12476718() {
        Path f1 = Paths.get("src/test/resources/ERR12476718/122-2_1.fq.gz");
        Path f2 = Paths.get("src/test/resources/ERR12476718/122-2_2.fq.gz");

        try {
            ValidatorWrapper vw = new ValidatorWrapper(
                    Arrays.asList(f1.toFile(), f2.toFile()), FileFormat.FASTQ, 100_000);
            vw.run();
            fail();
        } catch (ReadsValidationException e) {
            e.printStackTrace();
        }
    }
}
