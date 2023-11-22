package uk.ac.ebi.ena.readtools;

public class Read {
    protected final String name;
    protected final String bases;
    protected final String qualityScores;

    public Read(String name, String bases, String qualityScores) {
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
