package uk.ac.ebi.ena.readtools.refactored.provider;

import htsjdk.samtools.*;
import uk.ac.ebi.ena.readtools.cram.ref.ENAReferenceSource;
import uk.ac.ebi.ena.readtools.refactored.read.SamRead;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.Iterator;
import java.util.NoSuchElementException;

public class SamReadsProvider implements ReadsProvider<SamRead> {
    private final SamReader reader;
    private final Iterator<SAMRecord> samIterator;

    public SamReadsProvider(File samFile) {
        if (isCram(samFile)) {
            reader = SamReaderFactory.makeDefault()
                    .referenceSource(new ENAReferenceSource())
                    .open(samFile);
        } else {
            reader = SamReaderFactory.makeDefault()
                    .open(samFile);
        }
        samIterator = reader.iterator();
    }

    @Override
    public Iterator<SamRead> iterator() {
        return new Iterator<SamRead>() {
            @Override
            public boolean hasNext() {
                return samIterator.hasNext();
            }

            @Override
            public SamRead next() {
                if (!hasNext()) {
                    throw new NoSuchElementException();
                }

                SAMRecord samRecord = samIterator.next();
                return new SamRead(samRecord);
            }
        };
    }

    @Override
    public void close() throws Exception {
        if (reader != null) {
            reader.close();
        }
    }

    public static boolean isCram(File file) {
        // Magic number for CRAM files
        final byte[] cramMagicNumber = new byte[]{'C', 'R', 'A', 'M'};
        byte[] fileMagicNumber = new byte[4];

        try (FileInputStream fis = new FileInputStream(file)) {
            if (fis.read(fileMagicNumber) != cramMagicNumber.length) {
                return false; // File is shorter than the length of the magic number
            }
        } catch (IOException e) {
            e.printStackTrace();
            return false;
        }

        // Compare the file's magic number with the CRAM magic number
        for (int i = 0; i < cramMagicNumber.length; i++) {
            if (fileMagicNumber[i] != cramMagicNumber[i]) {
                return false; // Mismatch found
            }
        }

        return true; // Magic numbers match, it's a CRAM file
    }
}
