package uk.ac.ebi.ena.readtools.refactored.read;

import htsjdk.samtools.SAMRecord;
import uk.ac.ebi.ena.readtools.refactored.read.IRead;

public class SamRead implements IRead {
    private final SAMRecord samRecord;

    public SamRead(SAMRecord samRecord) {
        this.samRecord = samRecord;
    }

    public String getName() {
        return samRecord.getReadName();
    }

    public String getBases() {
        return samRecord.getReadString();
    }

    public String getQualityScores() {
        return samRecord.getBaseQualityString();
    }

    public boolean isReadLevelQuality() {
        return samRecord.getBaseQualityString().equals("*");
    }

    public boolean hasQualityControlFlag() {
        int flags = samRecord.getFlags();
        return (flags & 0x200) != 0;
    }

    public SAMRecord getSamRecord() {
        return samRecord;
    }
}
