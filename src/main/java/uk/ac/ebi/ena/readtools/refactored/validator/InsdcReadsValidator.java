package uk.ac.ebi.ena.readtools.refactored.validator;

import uk.ac.ebi.ena.readtools.refactored.provider.ReadsProvider;
import uk.ac.ebi.ena.readtools.refactored.read.IRead;

import java.util.Iterator;
import java.util.Set;
import java.util.HashSet;

public class InsdcReadsValidator implements ReadsValidator<IRead> {
    public static final String IUPAC_CODES = "ACGTURYSWKMBDHVNacgturyswkmbdhv.-";
    private static final int MIN_QUALITY_SCORE = 30;

    public static String ERROR_NULL_READS = "Reads cannot be null";
    public static String ERROR_NO_READS = "Submitted files must contain a minimum of 1 sequence read";
    public static String ERROR_EMPTY_READ = "Submitted files must not contain any empty reads";
    public static String ERROR_NOT_IUPAC = "Reads must contain only valid IUPAC codes, " +
            "with no more than 50% of bases being non-AUTCG";
    public static String ERROR_QUALITY = "When submitted file contains base quality scores " +
            "then >= 50% of reads must have average quality >= 30";

    @Override
    public boolean validate(ReadsProvider<IRead> provider) throws ReadsValidationException {
        long readCount = 0;
        long highQualityReadCount = 0;

        if (provider == null) {
            throw new ReadsValidationException(ERROR_NULL_READS, readCount);
        }

        Iterator<IRead> iterator = provider.iterator();
        if (!iterator.hasNext()) {
            throw new ReadsValidationException(ERROR_NO_READS, readCount);
        }
        Set<Character> iupacSet = new HashSet<>();
        for (char c : IUPAC_CODES.toCharArray()) {
            iupacSet.add(c);
        }

        while (iterator.hasNext()) {
            IRead read = iterator.next();
            String bases = read.getBases();
            String qualityScores = read.getQualityScores();

            readCount++;

            if (bases == null || bases.isEmpty()) {
                throw new ReadsValidationException(ERROR_EMPTY_READ, readCount);
            }

            int nonIUPACCount = 0;
            for (char base : bases.toCharArray()) {
                if (!iupacSet.contains(base)) {
                    nonIUPACCount++;
                }
            }

            if (((double) nonIUPACCount / bases.length()) > 0.5) {
                throw new ReadsValidationException(ERROR_NOT_IUPAC, readCount);
            }

            if (!qualityScores.isEmpty()) {
                int totalQuality = 0;
                for (char q : qualityScores.toCharArray()) {
                    totalQuality += q - '!'; // Phred+33 0 at !
                }
                if ((double) totalQuality / qualityScores.length() <= MIN_QUALITY_SCORE) {
                    highQualityReadCount++;
                }
            }

            if ((double) highQualityReadCount / readCount < 0.5) {
                throw new ReadsValidationException(ERROR_QUALITY, readCount);
            }
        }

        return true;
    }
}