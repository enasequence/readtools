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
package uk.ac.ebi.ena.readtools.sampler;

import java.io.File;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;

import uk.ac.ebi.ena.readtools.sampler.intervals.QualInterval;
import uk.ac.ebi.ena.readtools.sampler.intervals.QualScoring;

public class 
Main
{

    public static void
    main( String[] args ) throws Exception
    {
        Params params = new Params();
        
        JCommander jc = new JCommander( params );
        
        try
        {
            jc.parse( args );
        } catch( Exception e )
        {
            jc.usage();
            System.exit( 1 );
            
        }
        
        
        if( null != params.split && 0 < params.split.length() )
        {
            SplitReader sr = new SplitReader( new File( params.file_name ), new File( params.split ) );
            if( sr.process() )
            {
                System.out.println( "spots processed: " + sr.getProcessedSpotsCount() );
            }
        } else
        {
            FastqFileReader ffr = new FastqFileReader( new File( params.file_name ) );
            if( ffr.process() )
            {
                QualInterval qi = ffr.printIntervals();
                System.out.printf( ">>Picked: %s\n", qi );
                QualScoring  qs = ffr.printScoring( qi.getOffset() );
                System.out.printf( ">>Picked: %s\n", qs );
                ffr.printDataRow3();
                ffr.printDataRow2( qi.getOffset() );
            }   
        }
        
        
    }
    
    
    static class
    Params
    {
        @Parameter( names = "--file", description = "file name to sample", required = true )
        String file_name;
        
        @Parameter( names = "-s", description = "split" )
        String split;
/*        
        
        @Parameter( names = "--id", description = "run id to sample" )
        String run_id;
*/        
    }
    
}
