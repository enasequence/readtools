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
package uk.ac.ebi.ena.readtools.refactored;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;
import static uk.ac.ebi.ena.readtools.refactored.validator.SamReadsValidator.ERROR_FLAG_512;
import static uk.ac.ebi.ena.readtools.refactored.validator.SamReadsValidator.ERROR_QUAL_FIELD;

import org.junit.Test;

import uk.ac.ebi.ena.readtools.refactored.provider.ReadsProvider;
import uk.ac.ebi.ena.readtools.refactored.read.SamRead;
import uk.ac.ebi.ena.readtools.refactored.validator.ReadsValidationException;
import uk.ac.ebi.ena.readtools.refactored.validator.ReadsValidator;
import uk.ac.ebi.ena.readtools.refactored.validator.SamReadsValidator;

public class SamReadsValidatorTest {
    @Test
    public void readLevelQualityAsterisk() throws ReadsValidationException {
        ReadsProvider mrp;

        try {
            mrp = new MockReadsProvider(
                    new MockReadsProvider.MockRead("r1", "AGTC", "@@@@", false));
            new SamReadsValidator().validate(mrp);
            fail();
        } catch (ReadsValidationException e) {
            assertEquals(ERROR_QUAL_FIELD, e.getErrorMessage());
        }

        mrp = new MockReadsProvider(
                new MockReadsProvider.MockRead("r1", "AGTC", "*", true));
        new SamReadsValidator().validate(mrp);
    }

    @Test
    public void readLevelQualityFlag512() {
        try {
            ReadsProvider mrp = new MockReadsProvider(
                    new MockReadsProvider.MockRead("r1", "AGTC", "@@@@"));
            ReadsValidator<SamRead> v = new SamReadsValidator();
            v.validate(mrp);
            fail();
        } catch (ReadsValidationException e) {
            assertEquals(ERROR_FLAG_512, e.getErrorMessage());
        }
    }
}
