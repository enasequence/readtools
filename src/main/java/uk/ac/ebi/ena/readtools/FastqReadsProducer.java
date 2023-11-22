package uk.ac.ebi.ena.readtools;

import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;

import java.io.File;
import java.util.Iterator;
import java.util.NoSuchElementException;

public class FastqReadsProducer implements ReadsProducer {
    private FastqReader reader;

    public FastqReadsProducer(File fastqFile) {
        this.reader = new FastqReader(fastqFile);
    }

    @Override
    public Iterator<Read> iterator() {
        return new Iterator<Read>() {
            private FastqRecord nextRecord = reader.hasNext() ? reader.next() : null;

            @Override
            public boolean hasNext() {
                return nextRecord != null;
            }

            @Override
            public Read next() {
                if (nextRecord == null) {
                    throw new NoSuchElementException();
                }

                Read currentRead = new Read(
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
