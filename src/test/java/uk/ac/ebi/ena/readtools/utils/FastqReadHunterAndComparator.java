package uk.ac.ebi.ena.readtools.utils;

import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import org.junit.Assert;
import org.junit.Test;

import java.io.File;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.function.BiConsumer;
import java.util.function.Function;

/**
 * Contains a few quick tests to find and compare reads between fastq files produced with two Cram2Fasq versions i.e.<br/>
 *  - old/low-level version with older version of htsjdk<br/>
 *  - new/high-level version with newer version of htsjdk
 */
public class FastqReadHunterAndComparator {

    private static final String SOURCE_UNPAIRED_PATH_STR = "C:\\Users\\mhaseeb\\Documents\\ebi\\ena\\data-files\\20357-peritoneum-synchr-metastasis-ready\\20357-peritoneum-synchr-metastasis-ready-rt.fastq";
    private static final String SOURCE_PAIRED01_PATH_STR = "C:\\Users\\mhaseeb\\Documents\\ebi\\ena\\data-files\\20357-peritoneum-synchr-metastasis-ready\\20357-peritoneum-synchr-metastasis-ready-rt_1.fastq";
    private static final String SOURCE_PAIRED02_PATH_STR = "C:\\Users\\mhaseeb\\Documents\\ebi\\ena\\data-files\\20357-peritoneum-synchr-metastasis-ready\\20357-peritoneum-synchr-metastasis-ready-rt_2.fastq";

    private static final String TARGET_UNPAIRED_PATH_STR = "C:\\Users\\mhaseeb\\Documents\\ebi\\ena\\data-files\\20357-peritoneum-synchr-metastasis-ready\\20357-peritoneum-synchr-metastasis-ready-rtm.fastq";
    private static final String TARGET_PAIRED01_PATH_STR = "C:\\Users\\mhaseeb\\Documents\\ebi\\ena\\data-files\\20357-peritoneum-synchr-metastasis-ready\\20357-peritoneum-synchr-metastasis-ready-rtm_1.fastq";
    private static final String TARGET_PAIRED02_PATH_STR = "C:\\Users\\mhaseeb\\Documents\\ebi\\ena\\data-files\\20357-peritoneum-synchr-metastasis-ready\\20357-peritoneum-synchr-metastasis-ready-rtm_2.fastq";

    @Test
    public void simpleReadTest() {
        FastqReader fastqReader = new FastqReader(new File(SOURCE_PAIRED01_PATH_STR));

        AtomicInteger count = new AtomicInteger(0);
        fastqReader.forEach(fastqRecord -> {
            count.incrementAndGet();

//            System.out.println(String.format("BaseQualityHeader : %s, BaseQualityString : %s, ReadName : %s, ReadString : %s, FastqString : %s\n",
//            fastqRecord.getBaseQualityHeader(),
//            fastqRecord.getBaseQualityString(),
//            fastqRecord.getReadName(),
//            fastqRecord.getReadString(),
//            fastqRecord.toFastQString()));
        });

        System.out.println("Records : " + count.intValue());

        fastqReader.close();
    }

    @Test
    public void recordCountMatchTest() {

        Function<List<String>, Integer> func = pathStrList -> {
            AtomicInteger totalCount = new AtomicInteger(0);
            pathStrList.forEach(sourcePathStr -> {
                FastqReader fastqReader = new FastqReader(new File(sourcePathStr));

                fastqReader.forEach(fastqRecord -> {
                    totalCount.incrementAndGet();
                });

                fastqReader.close();
            });

            return totalCount.intValue();
        };

        int srcCount = func.apply(Arrays.asList(SOURCE_UNPAIRED_PATH_STR, SOURCE_PAIRED01_PATH_STR, SOURCE_PAIRED02_PATH_STR));
        int tgtCount = func.apply(Arrays.asList(TARGET_UNPAIRED_PATH_STR, TARGET_PAIRED01_PATH_STR, TARGET_PAIRED02_PATH_STR));

        Assert.assertEquals(srcCount, tgtCount);
    }

    @Test
    public void pairedReadsHunt() {
        FastqReader srcReader1 = new FastqReader(new File(SOURCE_PAIRED01_PATH_STR));
        FastqReader srcReader2 = new FastqReader(new File(SOURCE_PAIRED02_PATH_STR));

        FastqReader tgtReader1 = new FastqReader(new File(TARGET_PAIRED01_PATH_STR));
        FastqReader tgtReader2 = new FastqReader(new File(TARGET_PAIRED02_PATH_STR));

        int matchCount = 0, mismatchCount = 0;

        int readCount = 0;
        for(FastqRecord srcRead1 : srcReader1) {
            ++readCount;

            FastqRecord srcRead2 = null;
            if (srcReader2.hasNext()) {
                srcRead2 = srcReader2.next();
            } else {
                System.out.println("src reader2 empty.");
                break;
            }

            FastqRecord tgtRead1 = null;
            if (tgtReader1.hasNext()) {
                tgtRead1 = tgtReader1.next();
            } else {
                System.out.println("tgt reader1 empty.");
                break;
            }

            FastqRecord tgtRead2 = null;
            if (tgtReader2.hasNext()) {
                tgtRead2 = tgtReader2.next();
            } else {
                System.out.println("tgt reader2 empty.");
                break;
            }

//            if (readCount % 100_000 != 0) {
//                continue;
//            }

            if (tgtRead1.getReadString().equals(srcRead1.getReadString()) && tgtRead1.getBaseQualityString().equals(srcRead1.getBaseQualityString())) {
                if (tgtRead2.getReadString().equals(srcRead2.getReadString()) && tgtRead2.getBaseQualityString().equals(srcRead2.getBaseQualityString())) {
                    System.out.println(String.format("possible pair match. srcRead1 : %s, srcRead2 : %s, tgtRead1 : %s, tgtRead2 : %s",
                            srcRead1.getReadName(), srcRead2.getReadName(), tgtRead1.getReadName(), tgtRead2.getReadName()));

                    ++matchCount;
                } else {
                    ++mismatchCount;
                }
            } else {
                ++mismatchCount;
            }
        }
        srcReader1.close();
        srcReader2.close();

        tgtReader1.close();
        tgtReader2.close();

        System.out.println(String.format("matches : %d, mismatches : %d", matchCount, mismatchCount));
    }

    @Test
    public void readsPairingByOrder() {
        BiConsumer<String, String> pairPathsConsumer = (pathStr1, pathStr2) -> {

            FastqReader reader1 = new FastqReader(new File(pathStr1));
            FastqReader reader2 = new FastqReader(new File(pathStr2));

            int read1Num = 0;
            int read2Num = 0;

            for(FastqRecord read1 : reader1) {
                ++read1Num;

                FastqRecord read2 = null;
                if (reader2.hasNext()) {
                    read2 = reader2.next();
                    ++read2Num;
                } else {
                    System.out.println("reader2 empty.");
                    break;
                }

                if (read1Num % 10_000 != 0) {
                    continue;
                }

                System.out.println(String.format("possible pair. read1Num : %d, read2Num : %d, read1 : %s, read2 : %s",
                        read1Num, read2Num, read1.getReadName(), read2.getReadName()));
            }
            reader1.close();
            reader2.close();

            System.out.println(String.format("read1Count : %d, read2Count : %d", read1Num, read2Num));
        };

        pairPathsConsumer.accept(SOURCE_PAIRED01_PATH_STR, SOURCE_PAIRED02_PATH_STR);
        //pairPathsConsumer.accept(TARGET_PAIRED01_PATH_STR, TARGET_PAIRED02_PATH_STR);
    }

    @Test
    public void pairedReadsNameMismatches() {
        BiConsumer<String, String> pairPathsConsumer = (pathStr1, pathStr2) -> {
            int matchedNamesCount = 0, mismatchedNamesCount = 0;

            FastqReader reader1 = new FastqReader(new File(pathStr1));
            FastqReader reader2 = new FastqReader(new File(pathStr2));

            int read1Num = 0;
            int read2Num = 0;

            for(FastqRecord read1 : reader1) {
                ++read1Num;

                FastqRecord read2 = null;
                if (reader2.hasNext()) {
                    read2 = reader2.next();
                    ++read2Num;
                } else {
                    System.out.println("reader2 empty.");
                    break;
                }

                String read1Name = read1.getReadName();
                read1Name = read1Name.substring(0, read1Name.length() - 1);

                String read2Name = read2.getReadName();
                read2Name = read2Name.substring(0, read2Name.length() - 1);

                if (read1Name.equals(read2Name)) {
                    ++matchedNamesCount;

                    if (mismatchedNamesCount > 0) {
                        System.out.println(String.format("match after a mismatch. read1Num : %d, read2Num : %d, read1 : %s, read2 : %s, matchedNames : %d, mismatchedNames : %d",
                                read1Num, read2Num, read1.getReadName(), read2.getReadName(), matchedNamesCount, mismatchedNamesCount));
                    }
                } else {
                    ++mismatchedNamesCount;

                    if (mismatchedNamesCount == 1) {
                        System.out.println(String.format("first mismatch. read1Num : %d, read2Num : %d, read1 : %s, read2 : %s, matchedNames : %d, mismatchedNames : %d",
                            read1Num, read2Num, read1.getReadName(), read2.getReadName(), matchedNamesCount, mismatchedNamesCount));
                    }
                }
            }
            reader1.close();
            reader2.close();

            System.out.println(String.format("matchedNames : %d, mismatchedNames : %d", matchedNamesCount, mismatchedNamesCount));
            System.out.println(String.format("read1Count : %d, read2Count : %d", read1Num, read2Num));
        };

        //pairPathsConsumer.accept(SOURCE_PAIRED01_PATH_STR, SOURCE_PAIRED02_PATH_STR);
        pairPathsConsumer.accept(TARGET_PAIRED01_PATH_STR, TARGET_PAIRED02_PATH_STR);
    }
}
