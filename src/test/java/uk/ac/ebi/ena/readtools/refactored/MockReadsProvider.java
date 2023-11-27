package uk.ac.ebi.ena.readtools.refactored;

import uk.ac.ebi.ena.readtools.refactored.provider.ReadsProvider;
import uk.ac.ebi.ena.readtools.refactored.read.SamRead;

import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

public class MockReadsProvider implements ReadsProvider<MockReadsProvider.MockRead> {
    private final List<MockRead> reads;

    public MockReadsProvider(MockRead... reads) {
        this.reads = Arrays.asList(reads);
    }

    @Override
    public Iterator<MockRead> iterator() {
        return reads.iterator();
    }

    @Override
    public void close() {}

    public static class MockRead extends SamRead {
        String name;
        String bases;
        String qualityScores;
        private final boolean qualityControlFlag;

        public MockRead(String name, String bases, String qualityScores) {
            this(name, bases, qualityScores, false);
        }

        public MockRead(String name, String bases, String qualityScores, boolean qualityControlFlag) {
            super(null);

            this.name = name;
            this.bases = bases;
            this.qualityScores = qualityScores;
            this.qualityControlFlag = qualityControlFlag;
        }

        @Override
        public String getName() {
            return name;
        }

        @Override
        public String getBases() {
            return bases;
        }

        @Override
        public String getQualityScores() {
            return qualityScores;
        }

        @Override
        public boolean hasQualityControlFlag() {return qualityControlFlag;}
    }
}
