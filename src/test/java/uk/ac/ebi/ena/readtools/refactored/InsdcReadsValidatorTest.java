package uk.ac.ebi.ena.readtools.refactored;

import org.junit.Test;
import uk.ac.ebi.ena.readtools.refactored.provider.ReadsProvider;
import uk.ac.ebi.ena.readtools.refactored.validator.InsdcReadsValidator;
import uk.ac.ebi.ena.readtools.refactored.validator.ReadsValidationException;
import uk.ac.ebi.ena.readtools.refactored.MockReadsProvider.MockRead;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;
import static uk.ac.ebi.ena.readtools.refactored.validator.InsdcReadsValidator.*;

public class InsdcReadsValidatorTest {
    @Test
    public void noReads() {
        try {
            ReadsProvider mrp = new MockReadsProvider();
            new InsdcReadsValidator().validate(mrp);
            fail();
        } catch (ReadsValidationException e) {
            assertEquals(ERROR_NO_READS, e.getErrorMessage());
        }
    }

    @Test
    public void emptyRead() {
        ReadsProvider mrp;

        try {
            mrp = new MockReadsProvider(
                    new MockRead("r1", null, "1234"));
            new InsdcReadsValidator().validate(mrp);
            fail();
        } catch (ReadsValidationException e) {
            assertEquals(ERROR_EMPTY_READ, e.getErrorMessage());
        }

        try {
            mrp = new MockReadsProvider(
                    new MockRead("r1", "", "1234"));
            new InsdcReadsValidator().validate(mrp);
            fail();
        } catch (ReadsValidationException e) {
            assertEquals(ERROR_EMPTY_READ, e.getErrorMessage());
        }
    }

    @Test
    public void notIupac() {
        try {
            ReadsProvider mrp = new MockReadsProvider(
                    new MockRead("r1", "AFFF", "1234"));
            new InsdcReadsValidator().validate(mrp);
            fail();
        } catch (ReadsValidationException e) {
            assertEquals(ERROR_NOT_IUPAC, e.getErrorMessage());
        }
    }

    @Test
    public void lowQuality() throws ReadsValidationException {
        ReadsProvider mrp = new MockReadsProvider(
                new MockRead("r1", "AGTC", ">>>>"));
        new InsdcReadsValidator().validate(mrp);

        try {
            mrp = new MockReadsProvider(
                    new MockRead("r1", "AGTC", "@@@@"));
            new InsdcReadsValidator().validate(mrp);
            fail();
        } catch (ReadsValidationException e) {
            assertEquals(ERROR_QUALITY, e.getErrorMessage());
        }
    }
}
