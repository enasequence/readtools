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
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import uk.ac.ebi.ena.readtools.loader.common.QualityNormalizer;
import uk.ac.ebi.ena.readtools.loader.common.QualityNormalizer.QualityNormaizationException;
import uk.ac.ebi.ena.readtools.loader.common.feeder.DataProducerException;
import uk.ac.ebi.ena.readtools.loader.common.feeder.ProducibleData;
import uk.ac.ebi.ena.readtools.loader.common.feeder.ProducibleDataChecker;

public class 
CompleteGenomicsRead extends CompleteGenomicsBase implements Serializable 
{
    public interface 
    Flags
    {
        int LeftHalfDnbNoMatches     = 0x1;
        int LeftHalfDnbMapOverflow   = 0x2;
        int RightHalfDnbNoMatches    = 0x4;
        int RightHalfDnbMapOverflow  = 0x8;
    }
    
    static final long serialVersionUID = 1L;
    
    @ProducibleData( method = "feedLine" )
    public ReadFlags flags; // Mapping characteristics of the DNBs, represented in bits within an integer. Individual flags described below. 
    public String    reads;// The base calls read from a single DNB, in an order specified in lib_DNB_[LIBRARY-ID].tsv. Base positions for which no information is available are denoted by "N" in the reads field.
    public String    scores; //Quality scores for reads. Each score is a Phred-like transformation of the error probability associated with a single base read. Base positions for which no information is available are assigned a score of 0.
    public long      index;
    public long      key;
    
    @ProducibleDataChecker
    public void
    checkFeed() throws QualityNormaizationException
    {
        normailzer.normalize( scores );
    }
    
    final QualityNormalizer normailzer = QualityNormalizer.SANGER;
    
    // default - no constraints
    public
    CompleteGenomicsRead( ReadHeader read_header )
    {
        super( read_header );
    }
    

    public Object 
    getKey()
    {
        return -1 == key ? null : key;
    }
    
    public long 
    getIndex()
    {
        return index;
    }

    
    // regexs
    static Pattern p_line = Pattern.compile( "^(\\d{1,2})\t([ACGTN]{70})\t([!-~]{70})$" );
    private static long read_index = 0;
    private static long read_key   = 0;
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //feeding methods. NOTE: always use "public" access modifier because of method accession speed issues!
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    public void 
    feedLine( InputStream is ) throws IOException, DataProducerException
    {
        String line = readLine( is );
        Matcher m   = p_line.matcher( line );
        if( !m.find() )
            throw new DataProducerException( line_no - 1, String.format( "Line [%s] does not match regexp", line ) );
        
        flags  = ReadFlags.forString( m.group( 1 ) );
        reads  = m.group( 2 );
        scores = m.group( 3 );
        
        int read_flags = getFlags( m.group( 1 ) );
        if( -1 == read_flags ) 
            throw new DataProducerException( line_no - 1, String.format( "Flags: %d", read_flags ) );
        
        if( ( read_flags & ( Flags.LeftHalfDnbMapOverflow | Flags.RightHalfDnbMapOverflow ) ) == ( Flags.LeftHalfDnbMapOverflow | Flags.RightHalfDnbMapOverflow )  
            || ( read_flags & ( Flags.LeftHalfDnbMapOverflow | Flags.RightHalfDnbNoMatches ) ) == ( Flags.LeftHalfDnbMapOverflow | Flags.RightHalfDnbNoMatches )   
            || ( read_flags & ( Flags.LeftHalfDnbNoMatches | Flags.RightHalfDnbMapOverflow ) ) == ( Flags.LeftHalfDnbNoMatches | Flags.RightHalfDnbMapOverflow )
            || ( read_flags & ( Flags.LeftHalfDnbNoMatches | Flags.RightHalfDnbNoMatches ) ) == ( Flags.LeftHalfDnbNoMatches | Flags.RightHalfDnbNoMatches ) ) 
            key = -1;
        else
            key = read_key ++;
        
        index = read_index ++;
    }
    
    
    int 
    getFlags( String value )
    {
        return 1 == value.length() ? value.charAt( 0 ) - '0' 
                                   : 2 == value.length() ? value.charAt( 0 ) * 10 + value.charAt( 1 ) - '0' 
                                                         : -1;
    }
    
    
    public String
    toString()
    {
        return new StringBuilder()
            .append( "index  = [" )
            .append( index )
            .append( "]\n" )
            .append( "key    = [" )
            .append( getKey() )
            .append( "]\n" )
            .append( "flags  = [" )
            .append( flags )
            .append( "]\n" )
            .append( "reads  = [" )
            .append( reads )
            .append( "], length = " )
            .append( reads.length() )
            .append( "\nscores = [" )
            .append( scores )
            .append( "], length = " )
            .append( scores.length() )
            .toString();
    }
}
