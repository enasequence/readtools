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
package uk.ac.ebi.ena.readtools.v2.validator;

import java.io.File;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import htsjdk.samtools.SAMException;

import uk.ac.ebi.ena.readtools.loader.fastq.FastqIterativeWriter;
import uk.ac.ebi.ena.readtools.loader.fastq.PairedRead;
import uk.ac.ebi.ena.readtools.v2.provider.ReadsProvider;
import uk.ac.ebi.ena.readtools.v2.read.FastqRead;
import uk.ac.ebi.ena.readtools.webin.cli.rawreads.BloomWrapper;
import uk.ac.ebi.ena.readtools.webin.cli.rawreads.DelegateIterator;
import uk.ac.ebi.ena.readtools.webin.cli.rawreads.RawReadsFile;

public class FastqReadsValidator extends ReadsValidator<FastqRead> {
    /*
    @ Each sequence identifier line starts with @
1    <instrument> Characters
    allowed:
    a-z, A-Z, 0-9 and
    underscore
2    Instrument ID
    <run number> Numerical Run number on instrument
3    <flowcell
    ID>
    Characters
    allowed:
    a-z, A-Z, 0-9
4    <lane> Numerical Lane number
5    <tile> Numerical Tile number
6    <x_pos> Numerical X coordinate of cluster
7    <y_pos> Numerical Y coordinate of cluster
SPACE HERE
8    <read> Numerical Read number. 1 can be single read or read 2 of pairedend
9    <is
    filtered>
    Y or N Y if the read is filtered, N otherwise
10    <control
    number>
    Numerical 0 when none of the control bits are on, otherwise it is
    an even number. See below.
11    <index
    sequence>
    ACTG Index sequence
    */
    //                                                          1        :  2   :    3       :   4  :  5   :   6   :  7          8 :  9 :  10         : 11
//    final static private Pattern p_casava_1_8_name = Pattern.compile( "^@([a-zA-Z0-9_-]+:[0-9]+:[a-zA-Z0-9]+:[0-9]+:[0-9]+:[0-9-]+:[0-9-]+) ([12]):[YN]:[0-9]*[02468]:[ACGTN]+$" );
    // relaxed regular expression
    final static Pattern pCasava18Name = Pattern.compile(
            "^(.+)( +|\\t+)([0-9]+):[YN]:[0-9]*[02468]($|:.*$)");

    final static private Pattern pQuals = Pattern.compile("^([!-~]*?)$"); //qualities
    private static final int PAIRING_THRESHOLD = 20;
    private ReadStyle readStyle = null; // Field to keep track of the read style
    private final long expectedSize = 100;
    private final long readLimit = 20;

    public FastqReadsValidator(long readCountLimit) {
        super(readCountLimit);
    }

    @Override
    public boolean validate(ReadsProvider<FastqRead> provider) throws ReadsValidationException {
        BloomWrapper duplicationsBloomWrapper = new BloomWrapper(expectedSize);
        Map<String, Integer> counts = new HashMap<>(100);

        long readCount = 0;
        try {
            for (FastqRead read : provider) {
                if (readCount >= readCountLimit) {
                    break;
                }

                readCount++;
                if (readCount == 1) {
                    determineReadStyle(read.getName()); // Determine style based on the first read
                }

                validateRead(read, readCount);

                duplicationsBloomWrapper.add(read.getName());
                if (duplicationsBloomWrapper.getPossibleDuplicates().contains(read.getName())) {
                    counts.put(read.getName(), counts.getOrDefault(read.getName(), 0) + 1);
                }
            }

            if (duplicationsBloomWrapper.hasPossibleDuplicates()) {
                for (Map.Entry<String, Integer> e : counts.entrySet()) {
                    if (e.getValue() > 0) {
                        throw new ReadsValidationException(
                                String.format("Multiple (%d) occurrences of read name \"%s\" at: %s\n",
                                        e.getValue(),
                                        e.getKey(),
                                        e.getValue()));
                    }
                }
            }

            return true;
        } catch (SAMException e) {
            throw new ReadsValidationException(e.getMessage(), readCount);
        }
    }

    private void determineReadStyle(String name) {
        Matcher casavaMatcher = pCasava18Name.matcher(name);
        readStyle = casavaMatcher.matches() ? ReadStyle.CASAVA18 : ReadStyle.FASTQ;
    }

    private void validateRead(FastqRead read, long readCount) throws ReadsValidationException {
        validateReadName(read.getName(), readCount);
        validateQualityScores(read.getName(), read.getQualityScores(), readCount);
    }

    private void validateReadName(String name, long readCount) throws ReadsValidationException {
        switch (readStyle) {
            case CASAVA18:
                Matcher casavaMatcher = pCasava18Name.matcher(name);
                if (!casavaMatcher.matches()) {
                    throw new ReadsValidationException("Invalid CASAVA 1.8 read name", readCount, name);
                }
                break;
            case FASTQ:
                // For FASTQ, no specific pattern, but you can add basic checks if needed
                if (name == null || name.trim().isEmpty()) {
                    throw new ReadsValidationException("Invalid FASTQ read name", readCount, name);
                }
                // Other basic validations for FASTQ can be added here if needed
                break;
            default:
                throw new ReadsValidationException("Read style not determined for read name", readCount, name);
        }
    }

    private void validateQualityScores(String readName, String qualityScores, long readCount) throws ReadsValidationException {
        Matcher matcher = pQuals.matcher(qualityScores);
        if (!matcher.matches()) {
            throw new ReadsValidationException("Invalid quality scores: " + qualityScores, readCount, readName);
        }
    }

    private Map<String, Set<String>> findAllDuplications(BloomWrapper duplications, int limit, RawReadsFile rf) {
        Map<String, Integer> counts = new HashMap<>(limit);
        Map<String, Set<String>> results = new LinkedHashMap<>(limit);

        String msg = "Verifying possible duplicates for file " + rf.getFilename();
//        log.info(msg);

        FastqIterativeWriter wrapper = new FastqIterativeWriter();
        wrapper.setFiles(new File[]{new File(rf.getFilename())});
        wrapper.setReadType(FastqIterativeWriter.READ_TYPE.SINGLE);
        wrapper.setReadLimit(readLimit);

        Iterator<String> read_name_iterator = new DelegateIterator<PairedRead, String>(wrapper.iterator()) {
            @Override
            public String convert(PairedRead obj) {
                return obj.forward.getName();
            }
        };

        long index = 1;

        while (read_name_iterator.hasNext()) {
            String read_name = read_name_iterator.next();
            if (duplications.getPossibleDuplicates().contains(read_name)) {
                counts.put(read_name, counts.getOrDefault(read_name, 0) + 1);
                Set<String> dlist = results.getOrDefault(read_name, new LinkedHashSet<>());
                dlist.add(rf.getFilename() + ", read " + index);
                results.put(read_name, dlist);
            }
            ++index;
        }

        return results.entrySet()
                .stream()
                //only read names occurring more than once are considered duplicates
                .filter(e -> counts.get(e.getKey()) > 1)
                .limit(limit)
                .collect(Collectors.toMap(e -> e.getKey(), e -> e.getValue(), (v1, v2) -> v1, LinkedHashMap::new));
    }

    public enum ReadStyle {FASTQ, CASAVA18}
}
