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

import htsjdk.samtools.SAMRecord;
import uk.ac.ebi.ena.readtools.loader.common.writer.Spot;

import java.io.Serializable;

/**
 * Holds raw unpaired read information.
 */
public class SamRecordWrapper extends Read {
    private final SAMRecord samRecord;

    public SamRecordWrapper(SAMRecord samRecord) {
        super(samRecord.getReadName(), samRecord.getReadString(), samRecord.getBaseQualityString());
        this.samRecord = samRecord;
    }

    public SAMRecord getSamRecord() {
        return samRecord;
    }

    public long getSizeBytes() {
        return name.getBytes().length + bases.getBytes().length + qualityScores.getBytes().length;
    }
}
