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
import java.io.FileInputStream;
import java.io.IOException;
import java.util.Iterator;
import java.util.NoSuchElementException;

import htsjdk.samtools.*;

import uk.ac.ebi.ena.readtools.cram.ref.ENAReferenceSource;
import uk.ac.ebi.ena.readtools.refactored.read.SamRead;

public class SamReadsProvider implements ReadsProvider<SamRead> {
    private final SamReader reader;
    private final Iterator<SAMRecord> samIterator;

    public SamReadsProvider(File samFile) {
        if (isCram(samFile)) {
            reader = SamReaderFactory.makeDefault()
                    .referenceSource(new ENAReferenceSource())
                    .open(samFile);
        } else {
            reader = SamReaderFactory.makeDefault()
                    .open(samFile);
        }
        samIterator = reader.iterator();
    }

    @Override
    public Iterator<SamRead> iterator() {
        return new Iterator<SamRead>() {
            @Override
            public boolean hasNext() {
                return samIterator.hasNext();
            }

            @Override
            public SamRead next() {
                if (!hasNext()) {
                    throw new NoSuchElementException();
                }

                SAMRecord samRecord = samIterator.next();
                return new SamRead(samRecord);
            }
        };
    }

    @Override
    public void close() throws Exception {
        if (reader != null) {
            reader.close();
        }
    }

    public static boolean isCram(File file) {
        // Magic number for CRAM files
        final byte[] cramMagicNumber = new byte[]{'C', 'R', 'A', 'M'};
        byte[] fileMagicNumber = new byte[4];

        try (FileInputStream fis = new FileInputStream(file)) {
            if (fis.read(fileMagicNumber) != cramMagicNumber.length) {
                return false; // File is shorter than the length of the magic number
            }
        } catch (IOException e) {
            e.printStackTrace();
            return false;
        }

        // Compare the file's magic number with the CRAM magic number
        for (int i = 0; i < cramMagicNumber.length; i++) {
            if (fileMagicNumber[i] != cramMagicNumber[i]) {
                return false; // Mismatch found
            }
        }

        return true; // Magic numbers match, it's a CRAM file
    }
}
