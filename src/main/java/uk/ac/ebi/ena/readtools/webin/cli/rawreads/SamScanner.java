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
package uk.ac.ebi.ena.readtools.webin.cli.rawreads;

import static uk.ac.ebi.ena.readtools.sam.Sam2Fastq.INCLUDE_NON_PF_READS;
import static uk.ac.ebi.ena.readtools.sam.Sam2Fastq.INCLUDE_NON_PRIMARY_ALIGNMENTS;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;

import htsjdk.samtools.*;
import htsjdk.samtools.cram.ref.ReferenceSource;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.SequenceUtil;

public class SamScanner {
    public void checkSamFile2(File samFile, File referenceFile) throws IOException {
        SamReaderFactory.setDefaultValidationStringency(ValidationStringency.LENIENT);

        ReferenceSource referenceSource = new ReferenceSource(referenceFile);
        final SamReader samReader = SamReaderFactory.makeDefault()
                .referenceSource(referenceSource)
                .open(samFile);
        SAMFileHeader samHeader = samReader.getFileHeader();

        int totalReadCount = 0;
        for (final SAMRecord currentRecord : samReader) {

            if (currentRecord.isSecondaryOrSupplementary() && !INCLUDE_NON_PRIMARY_ALIGNMENTS) {
                continue;
            }

            // Skip non-PF reads as necessary
            if (currentRecord.getReadFailsVendorQualityCheckFlag() && !INCLUDE_NON_PF_READS) {
                continue;
            }

            ++totalReadCount;

            String readName = currentRecord.getReadName();

            byte[] readBases = Arrays.copyOf(currentRecord.getReadBases(), currentRecord.getReadBases().length);
            byte[] baseQualities = currentRecord.getBaseQualityString().getBytes(StandardCharsets.UTF_8);

            int totalBaseCount = readBases.length;
            if (totalBaseCount == 0) {
                throw new IOException("The file contains an empty read.");
            }

            if (currentRecord.getReadNegativeStrandFlag()) {
                SequenceUtil.reverseComplement(readBases);
                SequenceUtil.reverseQualities(baseQualities);
            }
        }

        CloserUtil.close(samReader);
    }
    public void checkSamFile(File inputFile, File referenceFile) throws IOException {
        SamReaderFactory readerFactory = SamReaderFactory.makeDefault();

        if (inputFile.getName().endsWith(".cram")) {
            if (referenceFile == null || !referenceFile.exists() || !referenceFile.isFile()) {
                throw new IOException("CRAM file provided but reference file is invalid or not provided.");
            }

            // Create a ReferenceSource object from the reference file
            ReferenceSource referenceSource = new ReferenceSource(referenceFile);
            readerFactory = readerFactory.referenceSource(referenceSource);
        }

        try (SamReader reader = readerFactory.open(inputFile)) {
            SAMRecordIterator iterator = reader.iterator();

            if (!iterator.hasNext()) {
                throw new IOException("The file must contain a minimum of 1 sequence read.");
            }

            int validReadsCount = 0;
            int highQualityReadsCount = 0;
            int totalReadsCount = 0;
            while (iterator.hasNext()) {
                SAMRecord record = iterator.next();
                totalReadsCount++;

                String readString = record.getReadString();
                if (readString.equals("*")) {
                    throw new IOException("The file contains an empty read.");
                }

                long invalidBases = readString.chars()
                        .filter(c -> !"AUTCGautcg".contains(String.valueOf((char) c)))
                        .count();
                if (invalidBases > readString.length() / 2) {
                    throw new IOException("Read contains more than 50% invalid IUPAC codes.");
                }
                validReadsCount++;

                byte[] baseQualities = record.getBaseQualities();
                if (baseQualities.length > 0) {
                    double averageQuality = 0;
                    for (byte quality : baseQualities) {
                        averageQuality += quality;
                    }
                    averageQuality /= baseQualities.length;
                    if (averageQuality >= 30) {
                        highQualityReadsCount++;
                    }
                }
            }

            if (validReadsCount > 0 && (double) highQualityReadsCount / validReadsCount < 0.5) {
                throw new IOException("Less than 50% of reads have average quality >= 30.");
            }

            System.out.println("File passed all checks.");
        }
    }
}
