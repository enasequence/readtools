package uk.ac.ebi.ena.readtools.refactored.provider;

import htsjdk.samtools.SAMException;
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import uk.ac.ebi.ena.readtools.loader.common.writer.ReadWriterException;
import uk.ac.ebi.ena.readtools.refactored.read.FastqRead;
import uk.ac.ebi.ena.readtools.refactored.validator.ReadsValidationException;

import java.io.File;
import java.util.Iterator;
import java.util.NoSuchElementException;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class FastqReadsProvider implements ReadsProvider<FastqRead> {
    private FastqReader reader;

    public FastqReadsProvider(File fastqFile) throws ReadsValidationException {
        try {
            this.reader = new FastqReader(fastqFile);
        } catch (SAMException e) {
            throw new ReadsValidationException(e.getMessage(), 1);
        }
    }

    @Override
    public Iterator<FastqRead> iterator() {
        return new Iterator<FastqRead>() {
            private FastqRecord nextRecord = reader.hasNext() ? reader.next() : null;

            @Override
            public boolean hasNext() {
                return nextRecord != null;
            }

            @Override
            public FastqRead next() {
                if (nextRecord == null) {
                    throw new NoSuchElementException();
                }

                FastqRead currentRead = new FastqRead(
                        nextRecord.getReadName(), nextRecord.getReadString(), nextRecord.getBaseQualityString());
                nextRecord = reader.hasNext() ? reader.next() : null;

                return currentRead;
            }
        };
    }

    @Override
    public void close() throws Exception {
        if (reader != null) {
            reader.close();
            reader = null;
        }
    }
}
