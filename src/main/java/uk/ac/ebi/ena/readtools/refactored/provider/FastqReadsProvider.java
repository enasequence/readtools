package uk.ac.ebi.ena.readtools.refactored.provider;

import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import uk.ac.ebi.ena.readtools.refactored.read.FastqRead;
import uk.ac.ebi.ena.readtools.refactored.read.IRead;

import java.io.File;
import java.util.Iterator;
import java.util.NoSuchElementException;

public class FastqReadsProvider implements ReadsProvider<FastqRead> {
    private FastqReader reader;

    public FastqReadsProvider(File fastqFile) {
        this.reader = new FastqReader(fastqFile);
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
