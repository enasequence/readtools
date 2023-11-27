package uk.ac.ebi.ena.readtools.refactored.validator;

import uk.ac.ebi.ena.readtools.refactored.provider.ReadsProvider;
import uk.ac.ebi.ena.readtools.refactored.read.SamRead;

import java.util.Iterator;

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
