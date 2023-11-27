package uk.ac.ebi.ena.readtools.refactored.validator;

public class ReadsValidationException extends Exception {
    private final String errorMessage;
    private final long readIndex;

    public ReadsValidationException(String errorMessage, long readIndex) {
        super(errorMessage);
        this.errorMessage = errorMessage;
        this.readIndex = readIndex;
    }

    public String getErrorMessage() {
        return errorMessage;
    }

    public long getReadIndex() {
        return readIndex;
    }
}
