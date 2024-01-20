package uk.ac.ebi.ena.readtools.v2.provider;

import uk.ac.ebi.ena.readtools.v2.FileFormat;
import uk.ac.ebi.ena.readtools.v2.read.IRead;
import uk.ac.ebi.ena.readtools.v2.validator.ReadsValidationException;

import java.io.File;

public class ReadsProviderFactory {
    private final File file;
    private final FileFormat format;
    private final boolean normaliseFastqQualityScores;

    public ReadsProviderFactory(File file, FileFormat format) {
        this(file, format, false);
    }

    public ReadsProviderFactory(File file, FileFormat format, boolean normaliseFastqQualityScores) {
        this.file = file;
        this.format = format;
        this.normaliseFastqQualityScores = normaliseFastqQualityScores;
    }

    public ReadsProvider<? extends IRead>
    makeReadsProvider() throws ReadsValidationException {
        switch (format) {
            case FASTQ:
                return new FastqReadsProvider(file, normaliseFastqQualityScores);
            case BAM:
            case CRAM:
                return new SamReadsProvider(file);
            default:
                throw new ReadsValidationException("not implemented");
        }
    }
}