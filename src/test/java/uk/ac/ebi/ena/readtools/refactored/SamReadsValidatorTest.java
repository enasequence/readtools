package uk.ac.ebi.ena.readtools.refactored;

import org.junit.Test;
import uk.ac.ebi.ena.readtools.refactored.provider.ReadsProvider;
import uk.ac.ebi.ena.readtools.refactored.read.SamRead;
import uk.ac.ebi.ena.readtools.refactored.validator.ReadsValidationException;
import uk.ac.ebi.ena.readtools.refactored.validator.ReadsValidator;
import uk.ac.ebi.ena.readtools.refactored.validator.SamReadsValidator;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;
import static uk.ac.ebi.ena.readtools.refactored.validator.SamReadsValidator.ERROR_FLAG_512;
import static uk.ac.ebi.ena.readtools.refactored.validator.SamReadsValidator.ERROR_QUAL_FIELD;

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
