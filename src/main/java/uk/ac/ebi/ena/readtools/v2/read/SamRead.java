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
package uk.ac.ebi.ena.readtools.v2.read;

import htsjdk.samtools.SAMRecord;


public class SamRead implements IRead {
    private final SAMRecord samRecord;

    public SamRead(SAMRecord samRecord) {
        this.samRecord = samRecord;
    }

    public String getName() {
        return samRecord.getReadName();
    }

    public String getBases() {
        return samRecord.getReadString();
    }

    public String getQualityScores() {
        return samRecord.getBaseQualityString();
    }

    public boolean isReadLevelQuality() {
        return samRecord.getBaseQualityString().equals("*");
    }

    public boolean hasQualityControlFlag() {
        int flags = samRecord.getFlags();
        return (flags & 0x200) != 0;
    }

    public SAMRecord getSamRecord() {
        return samRecord;
    }
}
