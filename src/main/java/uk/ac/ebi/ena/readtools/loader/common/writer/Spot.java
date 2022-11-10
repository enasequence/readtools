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
package uk.ac.ebi.ena.readtools.loader.common.writer;

import java.util.Objects;

public interface Spot {
    long getBaseCount();
    long getSizeBytes();

    public static class SpotName {
        public String name;
        public String index;

        public SpotName() {
        }

        public SpotName(String name, String index) {
            this.name = name;
            this.index = index;
        }

        @Override
        public String toString() {
            return "SpotName{" +
                    "name='" + name + '\'' +
                    ", index='" + index + '\'' +
                    '}';
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;
            SpotName spotName = (SpotName) o;
            return Objects.equals(name, spotName.name) && Objects.equals(index, spotName.index);
        }

        @Override
        public int hashCode() {
            return Objects.hash(name, index);
        }
    }
}
