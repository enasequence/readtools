package uk.ac.ebi.ena.readtools.common.reads.normalizers.htsjdk;

import htsjdk.samtools.SAMUtils;
import uk.ac.ebi.ena.readtools.common.reads.QualityNormalizer;

/**
 * Converts printable qualities in Sanger fastq format to binary phred scores.
 */
public class StandardQualityNormalizer implements QualityNormalizer {

    @Override
    public void normalize(byte[] qualities) {
        SAMUtils.fastqToPhred(qualities);
    }
}
