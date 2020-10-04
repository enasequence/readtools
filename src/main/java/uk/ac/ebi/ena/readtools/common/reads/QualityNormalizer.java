package uk.ac.ebi.ena.readtools.common.reads;

public interface QualityNormalizer {

    /**
     * Normalize qualities in-place.
     *
     * @param qualities
     */
    void normalize(byte[] qualities);
}
