package uk.ac.ebi.ena.readtools.common.reads.normalizers.htsjdk;

import htsjdk.samtools.util.SolexaQualityConverter;
import uk.ac.ebi.ena.readtools.common.reads.QualityNormalizer;

public class IlluminaQualityNormalizer implements QualityNormalizer {

    @Override
    public void normalize(byte[] qualities) {
        SolexaQualityConverter.getSingleton().convertSolexa_1_3_QualityCharsToPhredBinary(qualities);
    }
}
