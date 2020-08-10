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
package uk.ac.ebi.ena.readtools.loader.cg.reads;

import java.io.IOException;
import java.io.InputStream;
import java.io.Serializable;
import java.util.Arrays;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import uk.ac.ebi.ena.readtools.loader.common.QualityNormalizer.QualityNormaizationException;
import uk.ac.ebi.ena.readtools.loader.common.feeder.DataFeederException;
import uk.ac.ebi.ena.readtools.loader.common.feeder.FeedableData;
import uk.ac.ebi.ena.readtools.loader.common.feeder.FeedableDataChecker;

public class 
CompleteGenomicsMap3 extends CompleteGenomicsBase implements Serializable 
{
    public interface 
    Flags
    {
        int LastDNBRecord = 0x1;
        int side          = 0x2;
        int strand        = 0x4;
    }
    
    static final long serialVersionUID = 1L;

    //flags chromosome offsetInChr gap1 gap2 gap3 weight mateRec
    
    @FeedableData( method = "feedLine" )
    public int flags; 
    public String chromosome;
    public int offset; 
    public Integer gaps[];
    public int weight;
    public int mateRec;
    public int batch_index;
    public long index;
    
    

    @FeedableDataChecker
    public void
    checkFeed() throws QualityNormaizationException
    {
        //normailzer.normalize( scores );
    }
    
    // default - no constraints
    public
    CompleteGenomicsMap3( ReadHeader map_header )
    {
        super( map_header );
    }
    
    
    public Object
    getKey()
    {
        return index;
    }

    //0\tchr18\t54911965\t-2\t0\t5\t(\t1        "^([0-7])\t([\\S]+)\t([\\d]+)\t(-?[\\d]+)\t(-?[\\d]+)\t(-?[\\d]+\t){3}([!-~])\t([\\d])$"
    static Pattern p_line    = Pattern.compile( "^([0-7])\t([\\S]+)\t([\\d]+)\t(-?[\\d]{1,2})\t(-?[\\d]{1,2})\t(-?[\\d]{1,2})\t([!-~])\t([\\d]+)$" );
    static long    read_index;
    static int     dnb_batch_index;
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //feeding methods. NOTE: always use "public" access modifier because of method accession speed issues!
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    public void 
    feedLine( InputStream is ) throws IOException, DataFeederException
    {
        String line = readLine( is );
        Matcher m   = p_line.matcher( line );
        if( !m.find() )
            throw new DataFeederException( line_no, String.format( "Line [%s] does not match regexp", line ) );


        flags      = m.group( 1 ).charAt( 0 ) - '0';
        chromosome = m.group( 2 );
        offset     = Integer.parseInt( m.group( 3 ) );
        
        gaps = new Integer[ 3 ];
        gaps[ 0 ] = Integer.parseInt( m.group( 4 ) );
        gaps[ 1 ] = Integer.parseInt( m.group( 5 ) );
        gaps[ 2 ] = Integer.parseInt( m.group( 6 ) );

        weight  = m.group( 7 ).charAt( 0 ) - '!';
        mateRec = Integer.parseInt( m.group( 8 ) );
        index       = read_index;
        batch_index = dnb_batch_index;
        
        
        read_index += flags & Flags.LastDNBRecord;
        dnb_batch_index += Flags.LastDNBRecord == ( flags & Flags.LastDNBRecord ) ? - dnb_batch_index : 1; 
    }
    

    
    public String
    toString()
    {
        return new StringBuilder()
            .append( "index       = [" )
            .append( index )
            .append( ", " )
            .append( batch_index )
            .append( "]\n" )
            .append( "flags       = [" )
            .append( flags )
            .append( "]\n" )
            .append( "chromosome  = [" )
            .append( chromosome )
            .append( "]\noffset      = [" )
            .append( offset )
            .append( "]\ngaps        = " )
            .append( Arrays.asList( gaps ) )
            .append( "\nweight      = [" )
            .append( weight )
            .append( "]\nmateRec     = [" )
            .append( mateRec )
            .append( "]" )
            .toString();
    }
}
