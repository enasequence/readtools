/*
* Copyright 2010-2020 EMBL - European Bioinformatics Institute
* Licensed under the Apache License, Version 2.0 (the "License"); you may not use this
* file except in compliance with the License. You may obtain a copy of the License at
* http://www.apache.org/licenses/LICENSE-2.0
* Unless required by applicable law or agreed to in writing, software distributed under the
* License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
* CONDITIONS OF ANY KIND, either express or implied. See the License for the
* specific language governing permissions and limitations under the License.
*/
package uk.ac.ebi.ena.readtools.loader.fastq;


import uk.ac.ebi.ena.readtools.loader.common.consumer.Spot;

import java.util.StringJoiner;

/**
 * Holds unpaired or paired read information.
 */
public class FastqSpot implements Spot {

    public final String name;

    public final DataSpot forward;
    public final DataSpot reverse;

    public FastqSpot(String name, DataSpot dataSpot) {
        if (name == null || name.isEmpty() ) {
            throw new IllegalArgumentException("Invalid read name.");
        }

        if (dataSpot == null) {
            throw new IllegalArgumentException("Read cannot be null.");
        }

        this.forward = dataSpot;
        this.reverse = null;
        this.name = name;
    }

    /**
     * Constructor for a paired read.
     *
     * @param name
     * @param forward
     * @param reverse
     */
    public FastqSpot(String name, DataSpot forward, DataSpot reverse) {
        if (name == null || name.isEmpty() ) {
            throw new IllegalArgumentException("Invalid read name.");
        }

        if (forward == null && reverse == null) {
            throw new IllegalArgumentException("Both reads cannot be null.");
        }

        this.forward = forward;
        this.reverse = reverse;

        this.name = name;
    }

    @Override
    public String getName() {
        return name;
    }

    public boolean isPaired() {
        return !(reverse == null || forward == null);
    }

    /**
     * Convenient method for getting the non-null unpaired read.
     *
     * @return the non-null unpaired read. Returns null if this is a paired read and both forward and reverse reads are present.
     */
    public DataSpot getUnpaired() {
        if (forward == null) {
            return reverse;
        } else if (reverse == null) {
            return forward;
        } else {
            return null;
        }
    }

    @Override
    public String toString() {
        return new StringJoiner(", ", FastqSpot.class.getSimpleName() + "[", "]")
                .add("name='" + name + "'")
                .add("first=" + forward)
                .add("second=" + reverse)
                .toString();
    }
}
