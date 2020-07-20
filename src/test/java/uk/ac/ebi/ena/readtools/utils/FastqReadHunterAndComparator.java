package uk.ac.ebi.ena.readtools.utils;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.StringUtil;
import org.junit.Assert;
import org.junit.Test;
import uk.ac.ebi.ena.readtools.cram.ref.ENAReferenceSource;

import java.io.File;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.TreeMap;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.function.BiConsumer;
import java.util.function.Consumer;
import java.util.function.Function;

/**
 * Contains a few quick tests to find and compare reads between fastq files produced with two Cram2Fasq versions i.e.<br/>
 *  - old low-level version with older version of htsjdk<br/>
 *  - new high-level version with newer version of htsjdk
 */
public class FastqReadHunterAndComparator {

    private static final String SAM_SOURCE = "C:\\Users\\mhaseeb\\Documents\\ebi\\ena\\data-files\\20357-peritoneum-synchr-metastasis-ready\\20357-peritoneum-synchr-metastasis-ready.cram";

    //prod
    private static final String PROD_SOURCE_UNPAIRED_PATH_STR = "C:\\Users\\mhaseeb\\Documents\\ebi\\ena\\data-files\\20357-peritoneum-synchr-metastasis-ready\\ERR1300854.fastq";
    private static final String PROD_SOURCE_PAIRED01_PATH_STR = "C:\\Users\\mhaseeb\\Documents\\ebi\\ena\\data-files\\20357-peritoneum-synchr-metastasis-ready\\ERR1300854_1.fastq";
    private static final String PROD_SOURCE_PAIRED02_PATH_STR = "C:\\Users\\mhaseeb\\Documents\\ebi\\ena\\data-files\\20357-peritoneum-synchr-metastasis-ready\\ERR1300854_2.fastq";

    //files created using the old low-level generator.
    private static final String SOURCE_UNPAIRED_PATH_STR = "C:\\Users\\mhaseeb\\Documents\\ebi\\ena\\data-files\\20357-peritoneum-synchr-metastasis-ready\\20357-peritoneum-synchr-metastasis-ready-rt.fastq";
    private static final String SOURCE_PAIRED01_PATH_STR = "C:\\Users\\mhaseeb\\Documents\\ebi\\ena\\data-files\\20357-peritoneum-synchr-metastasis-ready\\20357-peritoneum-synchr-metastasis-ready-rt_1.fastq";
    private static final String SOURCE_PAIRED02_PATH_STR = "C:\\Users\\mhaseeb\\Documents\\ebi\\ena\\data-files\\20357-peritoneum-synchr-metastasis-ready\\20357-peritoneum-synchr-metastasis-ready-rt_2.fastq";

    //files created using the new high-level generator.
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

    /**
     * Iterate both sets of file comparing reads and qualities of the mates in the pairs with each other.<br/>
     * The mates in the second pair are not searched. Comparison is based only on the order of the reads.
     */
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

            //match first mates in the source and target pairs
            if (tgtRead1.getReadString().equals(srcRead1.getReadString()) && tgtRead1.getBaseQualityString().equals(srcRead1.getBaseQualityString())
                    && //match second mates in the source and target pairs
                    tgtRead2.getReadString().equals(srcRead2.getReadString()) && tgtRead2.getBaseQualityString().equals(srcRead2.getBaseQualityString())) {

                System.out.println(String.format("possible pair match. srcRead1 : %s, srcRead2 : %s, tgtRead1 : %s, tgtRead2 : %s",
                        srcRead1.getReadName(), srcRead2.getReadName(), tgtRead1.getReadName(), tgtRead2.getReadName()));

                ++matchCount;
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

    /**
     * Compare read names of the pairs with each other to see if they match.
     */
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

                //remove the index from the end.
                String read1Name = read1.getReadName();
                read1Name = read1Name.substring(0, read1Name.length() - 1);

                //remove the index from the end.
                String read2Name = read2.getReadName();
                read2Name = read2Name.substring(0, read2Name.length() - 1);

                if (read1Name.equals(read2Name)) {
                    ++matchedNamesCount;

                    // it would make sense that read names would never match with each other after the first mismatch because
                    // of some problem that messed it up at the time of file generation but surprisingly, they continue
                    // to match afterwards.
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

    @Test
    public void lookupReadsInGeneratedData() {
        SamPairedReadsLookup samPairedReadsLookup = new SamPairedReadsLookup();
        samPairedReadsLookup.samPathStr = SAM_SOURCE;
//        samPairedReadsLookup.fastq1PathStr = PROD_SOURCE_PAIRED01_PATH_STR;
//        samPairedReadsLookup.fastq2PathStr = PROD_SOURCE_PAIRED02_PATH_STR;
//        samPairedReadsLookup.fastqUnpairedPathStr = PROD_SOURCE_UNPAIRED_PATH_STR;
//        samPairedReadsLookup.fastq1PathStr = SOURCE_PAIRED01_PATH_STR;
//        samPairedReadsLookup.fastq2PathStr = SOURCE_PAIRED02_PATH_STR;
//        samPairedReadsLookup.fastqUnpairedPathStr = SOURCE_UNPAIRED_PATH_STR;
        samPairedReadsLookup.fastq1PathStr = TARGET_PAIRED01_PATH_STR;
        samPairedReadsLookup.fastq2PathStr = TARGET_PAIRED02_PATH_STR;
        samPairedReadsLookup.fastqUnpairedPathStr = TARGET_UNPAIRED_PATH_STR;

        samPairedReadsLookup.run();
    }

    @Test
    public void runValidationStats() {
        ValidationStats vsProd = new ValidationStats(PROD_SOURCE_PAIRED01_PATH_STR, PROD_SOURCE_PAIRED02_PATH_STR, PROD_SOURCE_UNPAIRED_PATH_STR);
        vsProd.run();
        vsProd.printStats();

        ValidationStats vsOld = new ValidationStats(SOURCE_PAIRED01_PATH_STR, SOURCE_PAIRED02_PATH_STR, SOURCE_UNPAIRED_PATH_STR);
        vsOld.run();
        vsOld.printStats();

        ValidationStats vsNew = new ValidationStats(TARGET_PAIRED01_PATH_STR, TARGET_PAIRED02_PATH_STR, TARGET_UNPAIRED_PATH_STR);
        vsNew.run();
        vsNew.printStats();
    }
}

class SamPairedReadsLookup {
    public static final int SEARCH_GROUP_SIZE = 10_000;

    public static boolean INCLUDE_NON_PRIMARY_ALIGNMENTS = false;
    public static boolean INCLUDE_NON_PF_READS = false;

    public String samPathStr;
    public String fastq1PathStr, fastq2PathStr, fastqUnpairedPathStr;

    int samPairCount = 0, searchAttemptCount = 0, matchFoundCount = 0, totalExactMatchCount = 0, highestMatchesInASearch = 0;
    int firstMateMatchCount = 0, secondMateMatchCount = 0;

    Random random = new Random();

    public void run() {
        SamReader samReader = SamReaderFactory.makeDefault().referenceSource(new ENAReferenceSource()).open(Paths.get(samPathStr));

        List<SAMRecord[]> pairedReadsSearchGroup = new ArrayList<>(SEARCH_GROUP_SIZE);

        final Map<String, SAMRecord> firstSeenMates = new TreeMap<>();
        for (final SAMRecord currentRecord : samReader) {

            if (currentRecord.isSecondaryOrSupplementary() && !INCLUDE_NON_PRIMARY_ALIGNMENTS) continue;
            if (currentRecord.getReadFailsVendorQualityCheckFlag() && !INCLUDE_NON_PF_READS) continue;

            final String currentReadName = currentRecord.getReadName();

            final SAMRecord firstRecord = firstSeenMates.remove(currentReadName);
            if (firstRecord == null) {
                firstSeenMates.put(currentReadName, currentRecord);
                continue;
            }

            ++samPairCount;

            SAMRecord firstMate = firstRecord.getFirstOfPairFlag() ? firstRecord : currentRecord;
            SAMRecord secondMate = firstRecord.getFirstOfPairFlag() ? currentRecord : firstRecord;

            pairedReadsSearchGroup.add(new SAMRecord[]{firstMate, secondMate});

            if (pairedReadsSearchGroup.size() == SEARCH_GROUP_SIZE) {
                SAMRecord[] pair = pairedReadsSearchGroup.get(random.nextInt(pairedReadsSearchGroup.size()));
                searchPairInFastqFile(pair[0], pair[1]);
                pairedReadsSearchGroup.clear();
            }
        }

        SAMRecord[] pair = pairedReadsSearchGroup.get(random.nextInt(pairedReadsSearchGroup.size()));
        searchPairInFastqFile(pair[0], pair[1]);

        CloserUtil.close(samReader);

        System.out.println(String.format("search complete. sampairs : %d, searchesattempted : %d, matchedpairs : %d, totalexactmatches : %d (%.2f matches/search), highestmatchesinasearch : %d, firstmatematches : %d, secondmatematches : %d",
                samPairCount, searchAttemptCount, matchFoundCount, totalExactMatchCount, (double)totalExactMatchCount / (double)searchAttemptCount, highestMatchesInASearch, firstMateMatchCount, secondMateMatchCount));
    }

    private void searchPairInFastqFile(SAMRecord firstMate, SAMRecord secondMate) {
        System.out.println(String.format("attempting to match. sampairs : %d, searchesattempted : %d, matchedpairs : %d, totalexactmatches : %d (%.2f matches/search), highestmatchesinasearch : %d, firstmatematches : %d, secondmatematches : %d",
                samPairCount, searchAttemptCount, matchFoundCount, totalExactMatchCount, (double)totalExactMatchCount / (double)searchAttemptCount, highestMatchesInASearch, firstMateMatchCount, secondMateMatchCount));

        String mate1SearchRead = getReadString(firstMate);
        String mate1SearchBaseQualities = getBaseQualitiesString(firstMate);

        String mate2SearchRead = getReadString(secondMate);
        String mate2SearchBaseQualities = getBaseQualitiesString(secondMate);

        int exactMatchCount = 0;

        FastqReader reader1 = new FastqReader(new File(fastq1PathStr));
        FastqReader reader2 = new FastqReader(new File(fastq2PathStr));
        while (reader1.hasNext() && reader2.hasNext()) {
            FastqRecord read1 = reader1.next(), read2 = reader2.next();

            boolean firstMateMatch = false;
            if (mate1SearchRead.equals(read1.getReadString()) && mate1SearchBaseQualities.equals(read1.getBaseQualityString())) {
                ++firstMateMatchCount;

                firstMateMatch = true;
            }

            if (mate2SearchRead.equals(read2.getReadString()) && mate2SearchBaseQualities.equals(read2.getBaseQualityString())) {
                ++secondMateMatchCount;

                if (firstMateMatch) {
                    ++exactMatchCount;
                    if (exactMatchCount == 1) {
                        ++matchFoundCount;

                        System.out.println(String.format("pair matched. sampairs : %d, searchesattempted : %d, matchedpairs : %d, totalexactmatches : %d (%.2f matches/search), highestmatchesinasearch : %d, firstmatematches : %d, secondmatematches : %d",
                                samPairCount, searchAttemptCount, matchFoundCount, totalExactMatchCount, (double)totalExactMatchCount / (double)searchAttemptCount, highestMatchesInASearch, firstMateMatchCount, secondMateMatchCount));
                    }
                }
            }
        }
        reader1.close();
        reader2.close();

        totalExactMatchCount += exactMatchCount;
        highestMatchesInASearch = Math.max(highestMatchesInASearch, exactMatchCount);

        System.out.println(String.format("exact matches. thissearch: %d, allsearches : %d, highestmatchesinasearch : %d", exactMatchCount, totalExactMatchCount, highestMatchesInASearch));

        ++searchAttemptCount;
    }

    private String getReadString(SAMRecord samRecord) {
        if (samRecord.getReadNegativeStrandFlag()) {
            return SequenceUtil.reverseComplement(samRecord.getReadString());
        } else {
            return samRecord.getReadString();
        }
    }

    private String getBaseQualitiesString(SAMRecord samRecord) {
        if (samRecord.getReadNegativeStrandFlag()) {
            return StringUtil.reverseString(samRecord.getBaseQualityString());
        } else {
            return samRecord.getBaseQualityString();
        }
    }
}

/**
 * - number of paired reads
 * - number of unpaired reads
 * - number of each base
 * - number of each quality score
 * - number or generated read names
 */
class ValidationStats {

    private final String pairedPathStr1, pairedPathStr2, unpairedPathStr;

    public int paired1RecCount = 0, paired2RecCount = 0, unpairedRecCount = 0;

    public Map<String, Long> paired1BaseOccurenceCountMap = new TreeMap<>();
    public Map<String, Long> paired1BaseQualityOccurenceCountMap = new TreeMap<>();

    public Map<String, Long> paired2BaseOccurenceCountMap = new TreeMap<>();
    public Map<String, Long> paired2BaseQualityOccurenceCountMap = new TreeMap<>();

    public Map<String, Long> unpairedBaseOccurenceCountMap = new TreeMap<>();
    public Map<String, Long> unpairedBaseQualityOccurenceCountMap = new TreeMap<>();

    public ValidationStats(String pairedPathStr1, String pairedPathStr2, String unpairedPathStr) {
        this.pairedPathStr1 = pairedPathStr1;
        this.pairedPathStr2 = pairedPathStr2;
        this.unpairedPathStr = unpairedPathStr;
    }

    public void run() {
        callForEachRecord(pairedPathStr1, this::findPaired1Stats);
        callForEachRecord(pairedPathStr2, this::findPaired2Stats);
        callForEachRecord(unpairedPathStr, this::findUnpairedStats);
    }

    private void callForEachRecord(String pathStr, Consumer<FastqRecord> consumer) {
        FastqReader reader = new FastqReader(Paths.get(pathStr).toFile());
        for (FastqRecord record : reader) {
            consumer.accept(record);
        }
        reader.close();
    }

    private void findPaired1Stats(FastqRecord record) {
        ++paired1RecCount;

        record.getReadString().chars().forEach(chr -> {
            updateMap(new String(new char[] {(char) chr}), paired1BaseOccurenceCountMap);
        });
        record.getBaseQualityString().chars().forEach(chr -> {
            updateMap(new String(new char[] {(char) chr}), paired1BaseQualityOccurenceCountMap);
        });
    }

    private void findPaired2Stats(FastqRecord record) {
        ++paired2RecCount;

        record.getReadString().chars().forEach(chr -> {
            updateMap(new String(new char[] {(char) chr}), paired2BaseOccurenceCountMap);
        });
        record.getBaseQualityString().chars().forEach(chr -> {
            updateMap(new String(new char[] {(char) chr}), paired2BaseQualityOccurenceCountMap);
        });
    }

    private void findUnpairedStats(FastqRecord record) {
        ++unpairedRecCount;

        record.getReadString().chars().forEach(chr -> {
            updateMap(new String(new char[] {(char) chr}), unpairedBaseOccurenceCountMap);
        });
        record.getBaseQualityString().chars().forEach(chr -> {
            updateMap(new String(new char[] {(char) chr}), unpairedBaseQualityOccurenceCountMap);
        });
    }

    private void updateMap(String chr, Map<String, Long> map) {
        Long count = map.get(chr);
        if (count == null) {
            count = 0L;
        }

        map.put(chr, ++count);
    }

    public void printStats() {
        System.out.println(String.format("paired1records : %d, paired2records : %d, totalpairedrecords : %d, unpairedrecords : %d, totalrecords : %d",
                paired1RecCount, paired2RecCount, paired1RecCount + paired2RecCount,
                unpairedRecCount, paired1RecCount + paired2RecCount + unpairedRecCount));

        //printBasesAndQualityStats();
        printTotalBasesAndQualityOccurrenceCount();
    }

    private void printBasesAndQualityStats() {
        System.out.println("paired1 bases:");
        printMap(paired1BaseOccurenceCountMap);
        System.out.println("paired1 quality scores:");
        printMap(paired1BaseQualityOccurenceCountMap);
    }

    private void printTotalBasesAndQualityOccurrenceCount() {
        Map<String, Long> totalMap = new TreeMap<>(paired1BaseOccurenceCountMap);
        updateCountFrom(paired2BaseOccurenceCountMap, totalMap);
        updateCountFrom(unpairedBaseOccurenceCountMap, totalMap);

        System.out.println("total bases counts:");
        printMap(totalMap);

        totalMap = new TreeMap<>(paired1BaseQualityOccurenceCountMap);
        updateCountFrom(paired2BaseQualityOccurenceCountMap, totalMap);
        updateCountFrom(unpairedBaseQualityOccurenceCountMap, totalMap);

        System.out.println("total bases qualities counts:");
        printMap(totalMap);
    }

    private void updateCountFrom(Map<String, Long> src, Map<String, Long> tgt) {
        src.entrySet().forEach(entrySet -> {
            Long count = tgt.get(entrySet.getKey());
            if (count == null) {
                count = 0L;
            }

            tgt.put(entrySet.getKey(), count + entrySet.getValue());
        });
    }

    private void printMap(Map<String, Long> map) {
        map.entrySet().stream().forEach(entrySet -> {
            System.out.println(String.format("\t%s: %d", entrySet.getKey(), entrySet.getValue()));
        });
    }
}
