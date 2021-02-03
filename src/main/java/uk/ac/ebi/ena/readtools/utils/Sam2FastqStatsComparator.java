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
package uk.ac.ebi.ena.readtools.utils;

import java.nio.file.Paths;
import java.util.Map;
import java.util.TreeMap;
import java.util.function.Consumer;

import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;

public class Sam2FastqStatsComparator {

    /**
     *
     * @param args - 6 args total. First 3 are paths to first set of fastq files. Last 3 are paths to the second set of fastq files.
     */
    public static void main(String[] args) {

        System.out.println("Collecting data from first set of files . . .");

        ValidationStats vs1 = new ValidationStats(args[0], args[1], args[2]);
        vs1.run();

        System.out.println("Collecting data from second set of files . . .");

        ValidationStats vs2 = new ValidationStats(args[3], args[4], args[5]);
        vs2.run();

        Map<Character, Long> m1 = vs1.getTotalBasesOccurrenceCounts();
        Map<Character, Long> m2 = vs2.getTotalBasesOccurrenceCounts();
        try {
            assertContentsMatch(m1, m2);
        } catch (Exception ex) {
            System.out.println("Error matching bases counts.\nfirst set :");
            printMap(m1);
            System.out.println("second set :");
            printMap(m2);

            throw ex;
        }

        m1 = vs1.getTotalBaseQualitiesOccurrenceCounts();
        m2 = vs2.getTotalBaseQualitiesOccurrenceCounts();
        try {
            assertContentsMatch(m1, m2);
        } catch (Exception ex) {
            System.out.println("Error matching base qualities counts.\nfirst set:");
            printMap(m1);
            System.out.println("second set :");
            printMap(m2);

            throw ex;
        }

        System.out.println("Bases and qualities match across all files.");
    }

    private static void assertContentsMatch(Map<Character, Long> map1, Map<Character, Long> map2) {
        if (map1.size() != map2.size()) {
            throw new RuntimeException("Sizes do not match.");
        }

        map1.entrySet().forEach(entrySet -> {
            if (!map2.containsKey(entrySet.getKey())) {
                throw new RuntimeException("Mapping not found for : " + entrySet.getKey());
            }

            if (!entrySet.getValue().equals(map2.get(entrySet.getKey()))) {
                throw new RuntimeException("Values do not match. value1 : " + entrySet.getValue() + ", value2 : " + map2.get(entrySet.getKey()));
            }
        });
    }

    private static void printMap(Map<? extends Object, ? extends Object> map) {
        map.entrySet().stream().forEach(entrySet -> {
            System.out.println(String.format("\t%s: %s", entrySet.getKey(), entrySet.getValue()));
        });
    }
}

class ValidationStats {

    private final String pairedPathStr1, pairedPathStr2, unpairedPathStr;

    public int paired1RecCount = 0, paired2RecCount = 0, unpairedRecCount = 0;
    public int paired1GenReadNameCount = 0, paired2GenReadNameCount = 0, unpairedGenReadNameCount = 0;

    public Map<Character, Long> paired1BaseOccurenceCountMap = new TreeMap<>();
    public Map<Character, Long> paired1BaseQualityOccurenceCountMap = new TreeMap<>();

    public Map<Character, Long> paired2BaseOccurenceCountMap = new TreeMap<>();
    public Map<Character, Long> paired2BaseQualityOccurenceCountMap = new TreeMap<>();

    public Map<Character, Long> unpairedBaseOccurenceCountMap = new TreeMap<>();
    public Map<Character, Long> unpairedBaseQualityOccurenceCountMap = new TreeMap<>();

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

    public Map<Character, Long> getTotalBasesOccurrenceCounts() {
        Map<Character, Long> totalMap = new TreeMap<>(paired1BaseOccurenceCountMap);
        updateCountFrom(paired2BaseOccurenceCountMap, totalMap);
        updateCountFrom(unpairedBaseOccurenceCountMap, totalMap);

        return totalMap;
    }

    public Map<Character, Long> getTotalBaseQualitiesOccurrenceCounts() {
        Map<Character, Long> totalMap = new TreeMap<>(paired1BaseQualityOccurenceCountMap);
        updateCountFrom(paired2BaseQualityOccurenceCountMap, totalMap);
        updateCountFrom(unpairedBaseQualityOccurenceCountMap, totalMap);

        return totalMap;
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

        if (isReadNameGenerated(record.getReadName())) {
            ++paired1GenReadNameCount;
        }

        record.getReadString().chars().forEach(chr -> {
            updateMap((char)chr, paired1BaseOccurenceCountMap);
        });
        record.getBaseQualityString().chars().forEach(chr -> {
            updateMap((char)chr, paired1BaseQualityOccurenceCountMap);
        });
    }

    private void findPaired2Stats(FastqRecord record) {
        ++paired2RecCount;

        if (isReadNameGenerated(record.getReadName())) {
            ++paired2GenReadNameCount;
        }

        record.getReadString().chars().forEach(chr -> {
            updateMap((char)chr, paired2BaseOccurenceCountMap);
        });
        record.getBaseQualityString().chars().forEach(chr -> {
            updateMap((char)chr, paired2BaseQualityOccurenceCountMap);
        });
    }

    private void findUnpairedStats(FastqRecord record) {
        ++unpairedRecCount;

        if (isReadNameGenerated(record.getReadName())) {
            ++unpairedGenReadNameCount;
        }

        record.getReadString().chars().forEach(chr -> {
            updateMap((char)chr, unpairedBaseOccurenceCountMap);
        });
        record.getBaseQualityString().chars().forEach(chr -> {
            updateMap((char)chr, unpairedBaseQualityOccurenceCountMap);
        });
    }

    private boolean isReadNameGenerated(String readName) {
        try {
            String withoutPrefix = readName.split(" ")[1];

            String toParse = withoutPrefix;
            if (withoutPrefix.length() > 2 && withoutPrefix.charAt(withoutPrefix.length() - 2) == '/') {
                String withoutPairSuffix = withoutPrefix.substring(0, withoutPrefix.length() - 2);
                toParse = withoutPairSuffix;
            }

            Integer.parseInt(toParse);

            return true;
        } catch (NumberFormatException ex) {
            return false;
        } catch (Exception ex) {
            throw new RuntimeException("error processing read name : " + readName, ex);
        }
    }

    private void updateMap(Character chr, Map<Character, Long> map) {
        Long count = map.get(chr);
        if (count == null) {
            count = 0L;
        }

        map.put(chr, ++count);
    }

    private void updateCountFrom(Map<Character, Long> src, Map<Character, Long> tgt) {
        src.entrySet().forEach(entrySet -> {
            Long count = tgt.get(entrySet.getKey());
            if (count == null) {
                count = 0L;
            }

            tgt.put(entrySet.getKey(), count + entrySet.getValue());
        });
    }

    public void printStats() {
        System.out.println(String.format("paired1records : %d, paired2records : %d, totalpairedrecords : %d, unpairedrecords : %d, totalrecords : %d",
                paired1RecCount, paired2RecCount, paired1RecCount + paired2RecCount,
                unpairedRecCount, paired1RecCount + paired2RecCount + unpairedRecCount));

        System.out.println(String.format("paired1genreadnames : %d, paired2genreadnames : %d, unpairedgenreadnames : %d, totalgenreadnames : %d",
                paired1GenReadNameCount, paired2GenReadNameCount, unpairedGenReadNameCount,
                paired1GenReadNameCount + paired2GenReadNameCount + unpairedGenReadNameCount));

        printBasesAndQualityStats();
        printTotalBasesAndQualityOccurrenceCount();
    }

    private void printBasesAndQualityStats() {
        System.out.println("paired1 bases:");
        printMap(paired1BaseOccurenceCountMap);
        System.out.println("paired1 quality scores:");
        printMap(paired1BaseQualityOccurenceCountMap);

        System.out.println("paired2 bases:");
        printMap(paired2BaseOccurenceCountMap);
        System.out.println("paired2 quality scores:");
        printMap(paired2BaseQualityOccurenceCountMap);

        System.out.println("unpaired bases:");
        printMap(unpairedBaseOccurenceCountMap);
        System.out.println("unpaired quality scores:");
        printMap(unpairedBaseQualityOccurenceCountMap);
    }

    private void printTotalBasesAndQualityOccurrenceCount() {
        System.out.println("total bases counts:");
        printMap(getTotalBasesOccurrenceCounts());

        System.out.println("total bases qualities counts:");
        printMap(getTotalBaseQualitiesOccurrenceCounts());
    }

    private void printMap(Map<Character, Long> map) {
        map.entrySet().stream().forEach(entrySet -> {
            System.out.println(String.format("\t%s: %d", entrySet.getKey(), entrySet.getValue()));
        });
    }
}
