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
package uk.ac.ebi.ena.readtools.v2.provider;

import static htsjdk.samtools.SAMUtils.phredToFastq;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.Iterator;
import java.util.NoSuchElementException;

import htsjdk.samtools.SAMException;
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.util.FastqQualityFormat;

import uk.ac.ebi.ena.readtools.common.reads.QualityNormalizer;
import uk.ac.ebi.ena.readtools.utils.Utils;
import uk.ac.ebi.ena.readtools.v2.read.FastqRead;
import uk.ac.ebi.ena.readtools.v2.validator.ReadsValidationException;

public class FastqReadsProvider implements ReadsProvider<FastqRead> {
    private FastqReader reader;
    private FastqQualityFormat qualityFormat = null;
    private QualityNormalizer qualityNormalizer;
    private final boolean normaliseQualityScores;

    public FastqReadsProvider(File fastqFile) throws ReadsValidationException {
        this(fastqFile, false);
    }

    public FastqReadsProvider(File fastqFile, boolean normaliseQualityScores) throws ReadsValidationException {
        this.normaliseQualityScores = normaliseQualityScores;
        try {
            this.reader = new FastqReader(fastqFile);
            if (this.normaliseQualityScores) {
                qualityFormat = Utils.detectFastqQualityFormat(fastqFile.getAbsolutePath(), null);
                qualityNormalizer = Utils.getQualityNormalizer(qualityFormat);
            }
        } catch (SAMException e) {
            throw new ReadsValidationException(e.getMessage());
        }
    }

    public FastqQualityFormat getQualityFormat() {
        return qualityFormat;
    }

    @Override
    public Iterator<FastqRead> iterator() {
        return new Iterator<FastqRead>() {
            private FastqRecord nextRecord = reader.hasNext() ? reader.next() : null;

            @Override
            public boolean hasNext() {
                return nextRecord != null;
            }

            @Override
            public FastqRead next() {
                if (nextRecord == null) {
                    throw new NoSuchElementException();
                }

                String normalisedQualityString = nextRecord.getBaseQualityString();
                if (normaliseQualityScores || qualityFormat == FastqQualityFormat.Standard) {
                    byte[] normalisedQualityBytes = normalisedQualityString.getBytes(StandardCharsets.UTF_8);
                    qualityNormalizer.normalize(normalisedQualityBytes);
                    normalisedQualityString = phredToFastq(normalisedQualityBytes);
                }

                FastqRead currentRead = new FastqRead(
                        nextRecord.getReadName(), nextRecord.getReadString(), normalisedQualityString);
                nextRecord = reader.hasNext() ? reader.next() : null;

                return currentRead;
            }
        };
    }

    @Override
    public void close() throws Exception {
        if (reader != null) {
            reader.close();
            reader = null;
        }
    }
}
