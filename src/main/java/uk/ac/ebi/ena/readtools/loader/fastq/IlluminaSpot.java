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


import uk.ac.ebi.ena.readtools.loader.common.consumer.DataConsumable;

public class
IlluminaSpot implements DataConsumable
{
    //Must be private. DO NOT USE!
    private 
    IlluminaSpot()
    {
        ;
    }
    
    
    public static final int FORWARD = 0;
    public static final int REVERSE = 1;

    //NCBI coords start from 0!
    public int[] read_start;
    public int[] read_length;
    public String[] read_name;
    //TODO byte arrays?
    public String bases = "";
    public String quals = "";
    public String name  = "";
    
    
    public static IlluminaSpot 
    initPaired()
    {
        IlluminaSpot result = new IlluminaSpot();
        result.read_start   = new int[] { -1, -1 };
        result.read_length  = new int[] { -1, -1 };
        result.read_name    = new String[ 2 ]; 
        return result;
    }
    
    
    public static IlluminaSpot 
    initSingle()
    {
        IlluminaSpot result = new IlluminaSpot();
        result.read_start   = new int[] { -1 };
        result.read_length  = new int[] { -1 };
        result.read_name    = new String[ 1 ];
        return result;
    }
    
    
    public String
    toString()
    {
        StringBuilder result = new StringBuilder()
                               .append( name )
                               .append( '\n' )
                               .append( read_start[ FORWARD ] )
                               .append( ":" )
                               .append( read_length[ FORWARD ] );
        
        if( 2 == read_start.length )
            result.append( ", " )
                  .append( read_start[ REVERSE ] )
                  .append( ":" )
                  .append( read_length[ REVERSE ] );
                   
        result.append( '\n' )
              .append( bases )
              .append( '\n' )
              .append( quals );
        
        return result.toString();
    }
}
