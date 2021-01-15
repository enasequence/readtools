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
package uk.ac.ebi.ena.readtools.loader.fastq;

import java.io.Serializable;

import uk.ac.ebi.ena.readtools.loader.common.consumer.Spot;

/**
 * Holds raw unpaired read information.
 */
public class DataSpot implements Serializable, Spot {

    static final long serialVersionUID = 1L;

    public String readIndex;

    public String name; // name for bases
    public String bases;// bases
    public String quals;

    @Override
    public String getName() {
        return name;
    }

    @Override
    public long getBaseCount() {
        return bases.length();
    }

    public String toString() {
        return new StringBuilder()
            .append( "base_name = [" )
            .append(name)
            .append( "]\n" )
            .append( "bases = [" )
            .append( bases )
            .append( "], length = " )
            .append( null == bases ? "null" : bases.length() )
            .append( "\nquals = [" )
            .append( quals )
            .append( "], length = " )
            .append( null == quals ? "null" : quals.length() )
            .toString();
    }
}
