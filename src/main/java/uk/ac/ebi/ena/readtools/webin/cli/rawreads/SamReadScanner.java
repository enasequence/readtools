package uk.ac.ebi.ena.readtools.webin.cli.rawreads;

import uk.ac.ebi.ena.readtools.loader.common.writer.ReadWriterException;
import uk.ac.ebi.ena.readtools.loader.fastq.Read;
import uk.ac.ebi.ena.webin.cli.validator.message.ValidationMessage;
import uk.ac.ebi.ena.webin.cli.validator.message.ValidationOrigin;
import uk.ac.ebi.ena.webin.cli.validator.message.ValidationResult;

public class SamReadScanner extends InsdcStandardCheckingScanner {
    private static final String BAM_STAR = "*";

    public SamReadScanner(int printFreq) {
        super(printFreq);
    }

    public void write(Read read) throws ReadWriterException {
        super.write(read);

        if(record.isSecondaryOrSupplementary() )
            continue;

        if(record.getDuplicateReadFlag())
            continue;

        if(record.getReadString().equals( BAM_STAR ) && record.getBaseQualityString().equals(BAM_STAR))
            continue;

        if( record.getReadBases().length != record.getBaseQualities().length )
        {
            ValidationResult readResult = result.create( new ValidationOrigin( "read number", read_no ) );
            readResult.add( ValidationMessage.error( "Mismatch between length of read bases and qualities" ) );
        }

        paired.compareAndSet( false, record.getReadPairedFlag() );
    }
}
