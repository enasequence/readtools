package uk.ac.ebi.ena.readtools.refactored.validator;

import uk.ac.ebi.ena.readtools.refactored.provider.ReadsProvider;
import uk.ac.ebi.ena.readtools.refactored.read.IRead;

public interface ReadsValidator<T extends IRead> {
    boolean validate(ReadsProvider<T> provider) throws ReadsValidationException;
}
