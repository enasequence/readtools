/*
* Copyright 2010-2021 EMBL - European Bioinformatics Institute
* Licensed under the Apache License, Version 2.0 (the "License"); you may not use this
* file except in compliance with the License. You may obtain a copy of the License at
* http://www.apache.org/licenses/LICENSE-2.0
* Unless required by applicable law or agreed to in writing, software distributed under the
* License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
* CONDITIONS OF ANY KIND, either express or implied. See the License for the
* specific language governing permissions and limitations under the License.
*/
package uk.ac.ebi.ena.readtools.v2;

import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

import uk.ac.ebi.ena.readtools.v2.provider.ReadsProvider;
import uk.ac.ebi.ena.readtools.v2.read.SamRead;

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
