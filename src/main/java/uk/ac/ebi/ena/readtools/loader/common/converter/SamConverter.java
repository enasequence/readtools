package uk.ac.ebi.ena.readtools.loader.common.converter;

import htsjdk.samtools.*;
import uk.ac.ebi.ena.readtools.cram.ref.ENAReferenceSource;
import uk.ac.ebi.ena.readtools.loader.common.writer.ReadWriter;
import uk.ac.ebi.ena.readtools.loader.fastq.Read;
import uk.ac.ebi.ena.readtools.loader.fastq.SamRecordWrapper;
import uk.ac.ebi.ena.webin.cli.validator.message.ValidationMessage;

import java.io.EOFException;
import java.io.IOException;
import java.io.InputStream;
import java.util.Arrays;

public class SamConverter extends AbstractReadConverter<Read> {

    private final boolean isCram;
    private SAMRecordIterator samRecordIterator;

    public SamConverter(boolean isCram, InputStream istream, ReadWriter<Read, ?> writer) {
        super(istream, writer);
        this.isCram = isCram;
    }

    public SamConverter(boolean isCram, InputStream istream, ReadWriter<Read, ?> writer, Long readLimit) {
        super(istream, writer, readLimit);
        this.isCram = isCram;
    }

    @Override
    protected void begin() {
        try {
            SamInputResource samInputResource = SamInputResource.of(istream);

            SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);
            SamReaderFactory factory = SamReaderFactory.make();
            factory.enable(SamReaderFactory.Option.DONT_MEMORY_MAP_INDEX);
            factory.validationStringency(ValidationStringency.SILENT);
            factory.samRecordFactory(DefaultSAMRecordFactory.getInstance());

            if (isCram) {
                ENAReferenceSource reference_source = new ENAReferenceSource();
                factory.referenceSource(reference_source);
            }
            SamReader samReader = factory.open(samInputResource);
            samRecordIterator = samReader.iterator();
        } catch (Exception ex) {
            throw new ConverterException(ex);
        }
    }

    @Override
    public Read getNextSpotFromInputStream(InputStream inputStream) throws IOException {
        if (samRecordIterator.hasNext()) {
            SAMRecord samRecord = samRecordIterator.next();
            return new SamRecordWrapper(samRecord);
        } else {
            throw new EOFException();
        }
    }
}
