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
    Integer index1 = null;
    Integer index2 = null;

    public PairedFastqWriter(File tmp_root, int spill_page_size, long spill_page_size_bytes, long spill_abandon_limit_bytes) {
        super(tmp_root, spill_page_size, spill_page_size_bytes, spill_abandon_limit_bytes);
    }

    @Override
    public String
    getKey(Read spot) {
        return spot.key;
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
        String readIndexStr = (null == spot.index) ? spot.defaultReadIndex : spot.index;

        int readIndex = Integer.parseInt(readIndexStr) - 1;

        if (null == index1) {
            index1 = readIndex;
        } else {
            if (null == index2 && readIndex != index1) {
                index2 = readIndex;
            }
        }

        int appliedIndex;
        if (readIndex == index1) {
            appliedIndex = 0;
        } else if (readIndex == index2) {
            appliedIndex = 1;
        } else {
            throw new RuntimeException("Unexpected read index " + spot);
        }

        if (null == list.get(appliedIndex)) {
            list.set(appliedIndex, spot);
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
