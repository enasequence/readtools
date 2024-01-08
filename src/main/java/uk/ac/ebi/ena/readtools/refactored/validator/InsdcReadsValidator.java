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
package uk.ac.ebi.ena.readtools.refactored.validator;

import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

import htsjdk.samtools.SAMException;

import uk.ac.ebi.ena.readtools.refactored.provider.ReadsProvider;
import uk.ac.ebi.ena.readtools.refactored.read.IRead;

public class InsdcReadsValidator extends ReadsValidator<IRead> {
    public static final String IUPAC_CODES = "ACGTURYSWKMBDHVNacgturyswkmbdhv.-";
    private static final int MIN_QUALITY_SCORE = 30;

    public static String ERROR_NULL_READS = "Reads cannot be null";
    public static String ERROR_NO_READS = "Submitted files must contain a minimum of 1 sequence read";
    public static String ERROR_EMPTY_READ = "Submitted files must not contain any empty reads";
    public static String ERROR_READ_NAME_LENGTH = "Read name length exceeds 256 characters";
    public static String ERROR_NOT_IUPAC = "Reads must contain only valid IUPAC codes, " +
            "with no more than 50% of bases being non-AUTCG";
    public static String ERROR_QUALITY = "When submitted file contains base quality scores " +
            "then >= 50% of reads must have average quality >= 30";
    public static String INVALID_FILE = "Invalid file structure";

    public InsdcReadsValidator(long readCountLimit) {
        super(readCountLimit);
    }

    @Override
    public boolean validate(ReadsProvider<IRead> provider) throws ReadsValidationException {
        long readCount = 0;
        long highQualityReadCount = 0;

        if (provider == null) {
            throw new ReadsValidationException(ERROR_NULL_READS, readCount);
        }

        Iterator<IRead> iterator;
        try {
            iterator = provider.iterator();
        } catch (SAMException e) {
            throw new ReadsValidationException(INVALID_FILE, readCount);
        }

        if (!iterator.hasNext()) {
            throw new ReadsValidationException(ERROR_NO_READS, readCount);
        }
        Set<Character> iupacSet = new HashSet<>();
        for (char c : IUPAC_CODES.toCharArray()) {
            iupacSet.add(c);
        }

        while (iterator.hasNext()) {
            if (readCount >= readCountLimit) {
                break;
            }

            IRead read = iterator.next();
            String bases = read.getBases();
            String qualityScores = read.getQualityScores();

            readCount++;

            if (bases == null || bases.isEmpty()) {
                throw new ReadsValidationException(ERROR_EMPTY_READ, readCount);
            }

            if (read.getName().trim().length() > 256) {
                throw new ReadsValidationException(ERROR_READ_NAME_LENGTH, readCount);
            }

            int nonIUPACCount = 0;
            for (char base : bases.toCharArray()) {
                if (!iupacSet.contains(base)) {
                    nonIUPACCount++;
                }
            }

            if (((double) nonIUPACCount / bases.length()) > 0.5) {
                throw new ReadsValidationException(ERROR_NOT_IUPAC, readCount);
            }

            if (!qualityScores.isEmpty()) {
                int totalQuality = 0;
                for (char q : qualityScores.toCharArray()) {
                    totalQuality += q - '!'; // Phred+33 0 at !
                }
                if ((double) totalQuality / qualityScores.length() <= MIN_QUALITY_SCORE) {
                    highQualityReadCount++;
                }
            }

            if ((double) highQualityReadCount / readCount < 0.5) {
                throw new ReadsValidationException(ERROR_QUALITY, readCount);
            }
        }

        return true;
    }
}
