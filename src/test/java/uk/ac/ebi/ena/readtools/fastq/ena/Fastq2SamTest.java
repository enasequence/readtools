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

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.math.BigInteger;
import java.nio.file.Files;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import org.apache.commons.lang3.exception.ExceptionUtils;
import org.junit.Assert;
import org.junit.Ignore;
import org.junit.Test;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;

import uk.ac.ebi.ena.readtools.loader.common.FileCompression;
import uk.ac.ebi.ena.readtools.loader.common.InvalidBaseCharacterException;
import uk.ac.ebi.ena.readtools.loader.common.converter.ConverterException;
import uk.ac.ebi.ena.readtools.loader.common.writer.ReadWriterException;
import uk.ac.ebi.ena.readtools.loader.common.writer.ReadWriterMemoryLimitException;

/**
 * TODO The tests in here need to compare data between source fastq(s) and generated bam file as well.
 */
public class Fastq2SamTest {
    @Test
    public void singleFastqReadAndBaseCount()
            throws IOException, ConverterException, ReadWriterException, NoSuchAlgorithmException {

        Fastq2Sam.Params params = new Fastq2Sam.Params();
        params.tmp_root = System.getProperty("java.io.tmpdir");
        params.sample_name = "SM-001";
        params.data_file = Files.createTempFile(null, ".bam").toString();
        params.compression = FileCompression.NONE.name();
        params.files = Arrays.asList(
                new File(Fastq2SamTest.class.getClassLoader()
                        .getResource("fastq_spots_correct_paired_with_unpaired_1.txt").getFile())
                        .getAbsolutePath());

        Fastq2Sam fastq2Sam = new Fastq2Sam();
        fastq2Sam.create(params);

        Assert.assertTrue(new File(params.data_file).length() > 0);
        Assert.assertEquals(4, fastq2Sam.getTotalReadCount());
        Assert.assertEquals(404, fastq2Sam.getTotalBaseCount());
        Assert.assertEquals("9017afbcef3ff94d0281dc847aebb067", calculateFileMd5(new File(params.data_file)));
    }

    @Test
    public void pairedFastqReadAndBaseCount()
            throws IOException, ConverterException, ReadWriterException, NoSuchAlgorithmException {

        Fastq2Sam.Params params = new Fastq2Sam.Params();
        params.tmp_root = System.getProperty("java.io.tmpdir");
        params.sample_name = "SM-001";
        params.data_file = Files.createTempFile(null, ".bam").toString();
        params.compression = FileCompression.NONE.name();
        params.files = Arrays.asList(
                new File(Fastq2SamTest.class.getClassLoader()
                        .getResource("fastq_spots_correct_paired_with_unpaired_1.txt").getFile())
                        .getAbsolutePath(),
                new File(Fastq2SamTest.class.getClassLoader()
                        .getResource("fastq_spots_correct_paired_with_unpaired_2.txt").getFile())
                        .getAbsolutePath());

        Fastq2Sam fastq2Sam = new Fastq2Sam();
        fastq2Sam.create(params);

        Assert.assertTrue(new File(params.data_file).length() > 0);
        Assert.assertEquals(8, fastq2Sam.getTotalReadCount());
        Assert.assertEquals(808, fastq2Sam.getTotalBaseCount());
        Assert.assertEquals("01ce849441f1d3ac174ce6c2bb435849", calculateFileMd5(new File(params.data_file)));
    }

    @Ignore("Only run manually if needed.")
    @Test
    public void twoLargeFilesManualTest() throws IOException, ConverterException, ReadWriterException {
        Fastq2Sam.Params params = new Fastq2Sam.Params();
        params.tmp_root = System.getProperty("java.io.tmpdir");
        params.sample_name = "SM-001";
        params.data_file = Files.createTempFile(null, ".bam").toString();
        params.compression = FileCompression.GZ.name();
        params.files = Arrays.asList(
                new File(Fastq2SamTest.class.getClassLoader()
                        .getResource("F1.fastq.gz").getFile()).getAbsolutePath(),
                new File(Fastq2SamTest.class.getClassLoader()
                        .getResource("F2.fastq.gz").getFile()).getAbsolutePath());

        params.spill_page_size_bytes = 2L * 1024L * 1024L * 1024L;
        params.spill_abandon_limit_bytes = 10L * 1024L * 1024L * 1024L;

        Fastq2Sam fastq2Sam = new Fastq2Sam();
        fastq2Sam.create(params);

        Assert.assertTrue(new File(params.data_file).length() > 0);
    }

    @Test
    public void pairedFastqPairNumber3()
            throws IOException, ConverterException, ReadWriterException, NoSuchAlgorithmException {

        Fastq2Sam.Params params = new Fastq2Sam.Params();
        params.tmp_root = System.getProperty("java.io.tmpdir");
        params.sample_name = "SM-001";
        params.data_file = Files.createTempFile(null, ".bam").toString();
        params.compression = FileCompression.NONE.name();
        params.files = Arrays.asList(
                new File(Fastq2SamTest.class.getClassLoader()
                        .getResource("fastq_spots_read_num_eq3_1.fastq").getFile()).getAbsolutePath(),
                new File(Fastq2SamTest.class.getClassLoader()
                        .getResource("fastq_spots_read_num_eq3_2.fastq").getFile()).getAbsolutePath());

        Fastq2Sam fastq2Sam = new Fastq2Sam();
        fastq2Sam.create(params);

        Assert.assertTrue(new File(params.data_file).length() > 0);
        Assert.assertEquals(6, fastq2Sam.getTotalReadCount());
        Assert.assertEquals(906, fastq2Sam.getTotalBaseCount());
        Assert.assertEquals("01ce849441f1d3ac174ce6c2bb435849", calculateFileMd5(new File(params.data_file)));
    }

    @Test
    public void testMemoryPaging()
            throws IOException, ConverterException, ReadWriterException, NoSuchAlgorithmException {

        Fastq2Sam.Params params = new Fastq2Sam.Params();
        params.tmp_root = System.getProperty("java.io.tmpdir");
        params.sample_name = "SM-001";
        params.data_file = Files.createTempFile(null, ".bam").toString();
        params.compression = FileCompression.NONE.name();
        params.files = Arrays.asList(
                new File(Fastq2SamTest.class.getClassLoader()
                        .getResource("fastq_spots_correct_paired_with_unpaired_1.txt").getFile())
                        .getAbsolutePath(),
                new File(Fastq2SamTest.class.getClassLoader()
                        .getResource("fastq_spots_correct_paired_with_unpaired_2.txt").getFile())
                        .getAbsolutePath());

        params.spill_page_size_bytes = 400L;
        params.spill_abandon_limit_bytes = 10L * 1024L;

        Fastq2Sam fastq2Sam = new Fastq2Sam();
        fastq2Sam.create(params);

        Assert.assertTrue(new File(params.data_file).length() > 0);
        Assert.assertEquals(8, fastq2Sam.getTotalReadCount());
        Assert.assertEquals(808, fastq2Sam.getTotalBaseCount());
        Assert.assertEquals("01ce849441f1d3ac174ce6c2bb435849", calculateFileMd5(new File(params.data_file)));
    }

    @Test
    public void testMemoryPagingAbandonLimit() throws IOException, ConverterException, ReadWriterException {
        Fastq2Sam.Params params = new Fastq2Sam.Params();
        params.tmp_root = System.getProperty("java.io.tmpdir");
        params.sample_name = "SM-001";
        params.data_file = Files.createTempFile(null, ".bam").toString();
        params.compression = FileCompression.NONE.name();
        params.files = Arrays.asList(
                new File(Fastq2SamTest.class.getClassLoader()
                        .getResource("fastq_spots_correct_paired_with_unpaired_1a.txt").getFile())
                        .getAbsolutePath(),
                new File(Fastq2SamTest.class.getClassLoader()
                        .getResource("fastq_spots_correct_paired_with_unpaired_2a.txt").getFile())
                        .getAbsolutePath());

        params.spill_page_size_bytes = 200L;
        params.spill_abandon_limit_bytes = 201L;

        Fastq2Sam fastq2Sam = new Fastq2Sam();

        try {
            fastq2Sam.create(params);
            fail();
        } catch (ReadWriterMemoryLimitException e) {
            assertTrue(e.getMessage().contains("Temp memory limit"));
        }
    }

    @Test
    public void testCorrectPairedWithUnpaired() throws IOException, ConverterException, ReadWriterException {
        File inpFile1 = new File(Fastq2SamTest.class.getClassLoader()
                .getResource("fastq_spots_correct_paired_with_unpaired_1.txt").getFile());
        File inpFile2 = new File(Fastq2SamTest.class.getClassLoader()
                .getResource("fastq_spots_correct_paired_with_unpaired_2.txt").getFile());

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

        try(SamReader samReader = SamReaderFactory.makeDefault().open(outFile);) {
            Assert.assertEquals(params.sample_name, samReader.getFileHeader().getReadGroup("A").getSample());

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

        int expectedRecCount = fastqRecordMap.values().stream()
                .map(list -> list.size())
                .reduce((listSize1, listSize2) -> listSize1 + listSize2).get();

        Assert.assertEquals(expectedRecCount, recCount);
    }

    @Test
    public void testUracilFastqSingle() throws IOException, NoSuchAlgorithmException {
        Fastq2Sam.Params params = new Fastq2Sam.Params();
        params.tmp_root = System.getProperty("java.io.tmpdir");
        params.sample_name = "SM-001";
        params.compression = FileCompression.NONE.name();
        params.files = Arrays.asList(
                new File(Fastq2SamTest.class.getClassLoader().getResource("uracil-bases_1.fastq").getFile())
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
        Assert.assertTrue(Pattern.matches("[uU]{1,2}", e.getInvalidCharacters().stream()
                .map(String::valueOf).collect(Collectors.joining())));

        // with flag set to true.
        params.convertUracil = true;
        params.data_file = Files.createTempFile(null, ".bam").toString();

        //Since U or u is not an acceptable base character inside a BAM file, execution of following line without
        //any error would mean that Uracil base conversion did take place and there is no need for its verification.
        Fastq2Sam fastq2Sam = new Fastq2Sam();
        fastq2Sam.create(params);

        Assert.assertTrue(new File(params.data_file).length() > 0);
        Assert.assertEquals(2, fastq2Sam.getTotalReadCount());
        Assert.assertEquals(40, fastq2Sam.getTotalBaseCount());
        Assert.assertEquals("fdc0986ec2ef619fa382b1d06566ba73", calculateFileMd5(new File(params.data_file)));
    }

    @Test
    public void testUracilFastqPaired() throws IOException, NoSuchAlgorithmException {
        Fastq2Sam.Params params = new Fastq2Sam.Params();
        params.tmp_root = System.getProperty("java.io.tmpdir");
        params.sample_name = "SM-001";
        params.compression = FileCompression.NONE.name();
        params.files = Arrays.asList(
                new File(Fastq2SamTest.class.getClassLoader().getResource("uracil-bases_1.fastq").getFile())
                        .getAbsolutePath(),
                new File(Fastq2SamTest.class.getClassLoader().getResource("uracil-bases_2.fastq").getFile())
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
        Assert.assertTrue(Pattern.matches("[uU]{1,2}", e.getInvalidCharacters().stream()
                .map(String::valueOf).collect(Collectors.joining())));

        // with flag set to true.
        params.convertUracil = true;
        params.data_file = Files.createTempFile(null, ".bam").toString();

        //Since U or u is not an acceptable base character inside a BAM file, execution of following line without
        //any error would mean that Uracil base conversion did take place and there is no need for its verification.
        Fastq2Sam fastq2Sam = new Fastq2Sam();
        fastq2Sam.create(params);

        Assert.assertTrue(new File(params.data_file).length() > 0);
        Assert.assertEquals(4, fastq2Sam.getTotalReadCount());
        Assert.assertEquals(80, fastq2Sam.getTotalBaseCount());
        Assert.assertEquals("9991e990fad7c39578be55e0ef12d6ac", calculateFileMd5(new File(params.data_file)));
    }

    private Map<String, List<FastqRecord>> createFastqRecordMap(File file1, File file2) {
        Map<String, List<FastqRecord>> fastqRecordMap = new LinkedHashMap<>();

        try (FastqReader reader = new FastqReader(file1)) {
            for (FastqRecord rec : reader) {
                fastqRecordMap.put(rec.getReadName().substring(0, rec.getReadName().length() - 2), new ArrayList<>(Arrays.asList(rec)));
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
                    //The read data in both mates should not be similar or there is no point of this test otherwise.
                    Assert.assertNotEquals(rec.getReadString(), val.get(0).getReadString());
                    val.add(rec);
                }
            }
        }

        return fastqRecordMap;
    }

    public static String calculateFileMd5(File file) throws IOException, NoSuchAlgorithmException {
        MessageDigest digest = MessageDigest.getInstance("MD5");
        byte[] buf = new byte[4096];
        int read = 0;
        try (BufferedInputStream is = new BufferedInputStream(new FileInputStream(file))) {
            while ((read = is.read(buf)) > 0)
                digest.update(buf, 0, read);

            byte[] message_digest = digest.digest();
            BigInteger value = new BigInteger(1, message_digest);
            return String.format(String.format("%%0%dx", message_digest.length << 1), value);
        }
    }
}
