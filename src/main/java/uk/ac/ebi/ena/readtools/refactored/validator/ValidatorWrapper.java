package uk.ac.ebi.ena.readtools.refactored.validator;

import uk.ac.ebi.ena.readtools.refactored.FileFormat;
import uk.ac.ebi.ena.readtools.refactored.provider.FastqReadsProvider;
import uk.ac.ebi.ena.readtools.refactored.provider.ReadsProvider;
import uk.ac.ebi.ena.readtools.refactored.provider.SamReadsProvider;

import java.io.File;
import java.util.List;

public class ValidatorWrapper {
    protected final List<File> files;
    protected final FileFormat format;

    public ValidatorWrapper(List<File> files, FileFormat format) {
        this.files = files;
        this.format = format;
    }

    public static void validateFastq(File file) throws ReadsValidationException {
        try (ReadsProvider producer = new FastqReadsProvider(file)) {
            new InsdcReadsValidator().validate(producer);
        } catch (ReadsValidationException rve) {
            throw rve;
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }

    public static void validateSam(File file) throws ReadsValidationException {
        try (ReadsProvider producer = new SamReadsProvider(file)) {
            new InsdcReadsValidator().validate(producer);
            new SamReadsValidator().validate(producer);
        } catch (ReadsValidationException rve) {
            throw rve;
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }
}
