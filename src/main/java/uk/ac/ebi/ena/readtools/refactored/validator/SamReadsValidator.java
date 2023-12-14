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

import java.util.Iterator;

import uk.ac.ebi.ena.readtools.refactored.provider.ReadsProvider;
import uk.ac.ebi.ena.readtools.refactored.read.SamRead;

public class SamReadsValidator implements ReadsValidator<SamRead> {
    public static String ERROR_QUAL_FIELD =
            "Read-level quality reads using BAM and CRAM files must have QUAL field value ‘*’";
    public static String ERROR_FLAG_512 =
            "Reads that do not pass quality control filters must set the FLAG 512 (0x200)";

    @Override
    public boolean validate(ReadsProvider<SamRead> provider) throws ReadsValidationException {
        Iterator<SamRead> iterator = provider.iterator();

        long readCount = 0;
        while (iterator.hasNext()) {
            SamRead read = iterator.next();
            readCount++;

            // SAM format specific checks
            validateSamRead(read, readCount);
        }

        return true;
    }

    private void validateSamRead(SamRead read, long readIndex) throws ReadsValidationException {
        // Check for specific SAM file validations
        if (!read.getQualityScores().equals("*")) {
            throw new ReadsValidationException(ERROR_QUAL_FIELD, readIndex);
        }

        if (!read.hasQualityControlFlag()) {
            throw new ReadsValidationException(ERROR_FLAG_512, readIndex);
        }
    }
}
