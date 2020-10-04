package uk.ac.ebi.ena.readtools.common.reads.normalizers.htsjdk;

import htsjdk.samtools.util.SolexaQualityConverter;
import uk.ac.ebi.ena.readtools.common.reads.QualityNormalizer;

public class SolexaQualityNormalizer implements QualityNormalizer {

    @Override
    public void normalize(byte[] qualities) {
        SolexaQualityConverter.getSingleton().convertSolexaQualityCharsToPhredBinary(qualities);
    }
}
