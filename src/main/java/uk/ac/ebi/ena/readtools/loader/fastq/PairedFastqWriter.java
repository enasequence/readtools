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
package uk.ac.ebi.ena.readtools.loader.fastq;

import java.io.File;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import uk.ac.ebi.ena.readtools.loader.common.writer.AbstractPagedReadWriter;
import uk.ac.ebi.ena.readtools.loader.common.writer.ReadWriterException;


public class
PairedFastqWriter extends AbstractPagedReadWriter<Read, PairedRead> {
    // Provided readname structure is @{readkey}{separator:1(.|/|:|_)}{index:1(0:1:2)}
    static final Pattern SPLIT_REGEXP = Pattern.compile("^(.*)(?:[\\.|:|/|_])([0-9]+)$");
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
    */
    //    A00953:544:HMTFHDSX3:2:1101:6768:1
    //             1        :  2   :    3       :   4  :  5   :   6   :  7
    //    "^([a-zA-Z0-9_-]+:[0-9]+:[a-zA-Z0-9]+:[0-9]+:[0-9]+:[0-9-]+:[0-9-]+)$"
    static final Pattern CASAVA_LIKE_EXCLUDE_REGEXP = Pattern.compile("^([a-zA-Z0-9_-]+:[0-9]+:[a-zA-Z0-9_-]+:[0-9]+:[0-9]+:[0-9-]+:[0-9-]+)$");

    static public final int KEY = 1;
    static public final int INDEX = 2;

    Integer index1 = null, index2 = null;

    public PairedFastqWriter(File tmp_root, int spill_page_size, long spill_page_size_bytes, long spill_abandon_limit_bytes) {
        super(tmp_root, spill_page_size, spill_page_size_bytes, spill_abandon_limit_bytes);
    }

    public static String
    getReadKey(String readname) throws ReadWriterException {
        return getReadPart(readname, KEY);
    }

    public static String
    getPairNumber(String readname) throws ReadWriterException {
        return getReadPart(readname, INDEX);
    }

    private static String
    getReadPart(String readname, int group) throws ReadWriterException {
        Matcher casavaLikeMatcher = CASAVA_LIKE_EXCLUDE_REGEXP.matcher(readname);
        if (!casavaLikeMatcher.find()) {
            Matcher m = SPLIT_REGEXP.matcher(readname);
            if (m.find()) {
                return m.group(group);
            }
        }
        throw new ReadWriterException(String.format("Readname [%s] does not match regexp", readname));
    }

    @Override
    public String
    getKey(Read spot) {
        try {
            return getReadKey(spot.name);
        } catch (ReadWriterException de) {
            return spot.name;
        }
    }

    public List<Read>
    newListBucket() {
        List<Read> list = super.newListBucket();
        list.add(null);
        list.add(null);
        return list;
    }

    public void
    append(List<Read> list, Read spot) throws ReadWriterException {
        String readIndexStr;
        try {
            readIndexStr = getPairNumber(spot.name);
        } catch (ReadWriterException de) {
            readIndexStr = spot.defaultReadIndex;
        }
        int readIndex = Integer.parseInt(readIndexStr) - 1;

        /**
         * Here read pair number found in the file
         * will be mapped to 1st and 2nd items in the list
         */
        if (null == index1) {
            index1 = readIndex;
        } else {
            if (null == index2) {
                index2 = readIndex;
            }
        }
        if (readIndex != index1 && readIndex != index2) {
            throw new ReadWriterException(
                    "Unexpected read pair number: " + readIndex
                            + "; pair numbers " + index1 + " and " + index2 + " were found previously in the file.");
        }
        int mappedIndex = (readIndex == index1) ? 0 : 1;

        if (null == list.get(mappedIndex)) {
            list.set(mappedIndex, spot);
        } else {
            throw new RuntimeException("Got same spot twice: " + spot);
        }
    }

    @Override
    public PairedRead
    assemble(final String key, List<Read> list) throws ReadWriterException {
        PairedRead spot = list.size() == 1
                ? new PairedRead(key, list.get(0))
                : new PairedRead(key, list.get(0), list.get(1));

        return spot;
    }

    @Override
    public boolean
    isCollected(List<Read> list) {
        return null != list.get(0)
                && null != list.get(1);
    }

    @Override
    public PairedRead
    handleErrors(final String key, List<Read> list) throws ReadWriterException {
        return assemble(key, list);
    }
}
