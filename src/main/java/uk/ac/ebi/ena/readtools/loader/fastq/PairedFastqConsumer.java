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

import uk.ac.ebi.ena.readtools.loader.common.consumer.AbstractPagedDataConsumer;
import uk.ac.ebi.ena.readtools.loader.common.consumer.DataConsumerException;


public class
PairedFastqConsumer extends AbstractPagedDataConsumer<DataSpot, FastqSpot> {
    // Provided readname structure is @{readkey}{separator:1(.|/|:|_)}{index:1(0:1:2)}
    static final Pattern split_regexp = Pattern.compile("^(.*)(?:[\\.|:|/|_])([12])$");
    static public final int KEY = 1;
    static public final int INDEX = 2;

    public PairedFastqConsumer(File tmp_root, int spill_page_size, long spill_page_size_bytes) {
        super(tmp_root, spill_page_size, spill_page_size_bytes);
    }

    public static String
    getReadKey(String readname) throws DataConsumerException {
        return getReadPart(readname, KEY);
    }

    public static String
    getReadIndex(String readname) throws DataConsumerException {
        return getReadPart(readname, INDEX);
    }

    private static String
    getReadPart(String readname, int group) throws DataConsumerException {
        Matcher m = split_regexp.matcher(readname);
        if (m.find())
            return m.group(group);

        throw new DataConsumerException(String.format("Readname [%s] does not match regexp", readname));
    }

    @Override
    public String
    getKey(DataSpot spot) {
        try {
            return getReadKey(spot.name);
        } catch (DataConsumerException de) {
            return spot.name;
        }
    }

    public List<DataSpot>
    newListBucket() {
        List<DataSpot> list = super.newListBucket();
        list.add(null);
        list.add(null);
        return list;
    }

    public void
    append(List<DataSpot> list, DataSpot spot) throws DataConsumerException {
        String readIndexStr;
        try {
            readIndexStr = getReadIndex(spot.name);
        } catch (DataConsumerException de) {
            readIndexStr = spot.readIndex;
        }

        int readIndex = Integer.parseInt(readIndexStr) - 1;
        if (null == list.get(readIndex)) {
            list.set(readIndex, spot);
        } else {
            throw new RuntimeException("Got same spot twice: " + spot);
        }
    }

    @Override
    public FastqSpot
    assemble(final String key, List<DataSpot> list) throws DataConsumerException {
        FastqSpot spot = list.size() == 1
                ? new FastqSpot(key, list.get(0))
                : new FastqSpot(key, list.get(0), list.get(1));

        return spot;
    }

    @Override
    public boolean
    isCollected(List<DataSpot> list) {
        return null != list.get(0)
                && null != list.get(1);
    }

    @Override
    public FastqSpot
    handleErrors(final String key, List<DataSpot> list) throws DataConsumerException {
        return assemble(key, list);
    }
}
