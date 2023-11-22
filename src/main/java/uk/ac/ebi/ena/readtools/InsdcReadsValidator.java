package uk.ac.ebi.ena.readtools;

import java.util.Iterator;
import java.util.Set;
import java.util.HashSet;

public class InsdcReadsValidator {

    private static final String IUPAC_CODES = "ACGTURYSWKMBDHVN.-";
    private static final int MIN_QUALITY_SCORE = 30;

    public boolean validate(Iterable<Read> reads) {
        if (reads == null) {
            throw new IllegalArgumentException("Reads cannot be null");
        }

        Iterator<Read> iterator = reads.iterator();
        if (!iterator.hasNext()) {
            throw new IllegalArgumentException("Submitted files must contain a minimum of 1 sequence read");
        }

        int readCount = 0;
        int highQualityReadCount = 0;
        Set<Character> iupacSet = new HashSet<>();
        for (char c : IUPAC_CODES.toCharArray()) {
            iupacSet.add(c);
        }

        while (iterator.hasNext()) {
            Read read = iterator.next();
            String bases = read.getBases();
            String qualityScores = read.getQualityScores();

            if (bases.isEmpty()) {
                throw new IllegalArgumentException("Submitted files must not contain any empty reads");
            }

            int nonIUPACCount = 0;
            for (char base : bases.toCharArray()) {
                if (!iupacSet.contains(base)) {
                    nonIUPACCount++;
                }
            }

            if (((double) nonIUPACCount / bases.length()) > 0.5) {
                throw new IllegalArgumentException(
                        "Reads must contain only valid IUPAC codes, with no more than 50% of bases being non-AUTCG");
            }

            if (!qualityScores.isEmpty()) {
                int totalQuality = 0;
                for (char q : qualityScores.toCharArray()) {
                    totalQuality += q - '!';
                }
                if ((double) totalQuality / qualityScores.length() >= MIN_QUALITY_SCORE) {
                    highQualityReadCount++;
                }
            }

            readCount++;
        }

        if (!reads.iterator().next().getQualityScores().isEmpty() && (double) highQualityReadCount / readCount < 0.5) {
            throw new IllegalArgumentException(
                    "When submitted file contains base quality scores " +
                            "then >= 50% of reads must have average quality >= 30");
        }

        return true;
    }
}
