package uk.ac.ebi.ena.readtools.refactored.read;

import uk.ac.ebi.ena.readtools.refactored.read.IRead;

public class FastqRead implements IRead {
    private final String name;
    private final String bases;
    private final String qualityScores;

    public FastqRead(String name, String bases, String qualityScores) {
        this.name = name;
        this.bases = bases;
        this.qualityScores = qualityScores;
    }

    public String getName() {
        return name;
    }

    public String getBases() {
        return bases;
    }

    public String getQualityScores() {
        return qualityScores;
    }
}
