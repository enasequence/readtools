package uk.ac.ebi.ena.readtools.refactored;

import org.junit.Test;
import uk.ac.ebi.ena.readtools.refactored.provider.FastqReadsProvider;
import uk.ac.ebi.ena.readtools.refactored.provider.ReadsProvider;
import uk.ac.ebi.ena.readtools.refactored.validator.FastqReadsValidator;
import uk.ac.ebi.ena.readtools.refactored.validator.ReadsValidationException;
import uk.ac.ebi.ena.readtools.webin.cli.rawreads.FastqScannerTest;

import java.io.File;

public class FastqReadsValidatorTest {
    @Test
    public void t1() throws ReadsValidationException {
        ReadsProvider mrp = new FastqReadsProvider(
                new File(String.valueOf(FastqScannerTest.class.getClassLoader().getResource(
                        "rawreads/EP0_GTTCCTT_S1.txt.gz"))));
        new FastqReadsValidator().validate(mrp);
    }
}
