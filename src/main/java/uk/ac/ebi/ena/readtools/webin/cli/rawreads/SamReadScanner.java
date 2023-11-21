package uk.ac.ebi.ena.readtools.webin.cli.rawreads;

import htsjdk.samtools.SAMRecord;
import uk.ac.ebi.ena.readtools.loader.common.writer.ReadWriterException;
import uk.ac.ebi.ena.readtools.loader.fastq.Read;
import uk.ac.ebi.ena.readtools.loader.fastq.SamRecordWrapper;
import uk.ac.ebi.ena.webin.cli.validator.message.ValidationMessage;
import uk.ac.ebi.ena.webin.cli.validator.message.ValidationOrigin;
import uk.ac.ebi.ena.webin.cli.validator.message.ValidationResult;

public abstract class SamReadScanner extends InsdcStandardCheckingScanner {
    public SamReadScanner(int printFreq) {
        super(printFreq);
    }

    public void write(SamRecordWrapper read) throws ReadWriterException {
        if(read.getSamRecord().isSecondaryOrSupplementary()) {
            return;
        }
        if(read.getSamRecord().getDuplicateReadFlag()) {
            return;
        }

        super.write(read);
    }
}
