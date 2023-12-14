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
package uk.ac.ebi.ena.readtools.refactored.provider;

import java.io.File;
import java.util.Iterator;
import java.util.NoSuchElementException;

import htsjdk.samtools.SAMException;
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;

import uk.ac.ebi.ena.readtools.refactored.read.FastqRead;
import uk.ac.ebi.ena.readtools.refactored.validator.ReadsValidationException;

public class FastqReadsProvider implements ReadsProvider<FastqRead> {
    private FastqReader reader;

    public FastqReadsProvider(File fastqFile) throws ReadsValidationException {
        try {
            this.reader = new FastqReader(fastqFile);
        } catch (SAMException e) {
            throw new ReadsValidationException(e.getMessage(), 1);
        }
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

                FastqRead currentRead = new FastqRead(
                        nextRecord.getReadName(), nextRecord.getReadString(), nextRecord.getBaseQualityString());
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
