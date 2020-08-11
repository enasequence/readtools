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

import java.io.EOFException;
import java.io.IOException;
import java.io.InputStream;

public abstract class  
CompleteGenomicsBase
{
    
    CompleteGenomicsBase( ReadHeader header )
    {
        this.header = header;
    }
    
    public long line_no = 1;
    final public ReadHeader header;
    
    
    abstract Object getKey();
    
    // reads stream line till line separator
    protected String
    readLine( InputStream istream ) throws IOException
    {
        return readLine( istream, -1 );
    }

    
    // reads stream line till either the length or in case of len = -1 till line separator 
    protected String
    readLine( InputStream istream, 
              long len ) throws IOException
    {
            return readLine( istream, len, -1 );
    }
        

    // reads stream line till either the length or in case of len = -1 till line separator
    // here is unfair checking of stop_seq ( only 1 symbol is allowed!!! )
    protected String
    readLine( InputStream istream, 
              long        len,
              int         stop ) throws IOException
    {
        StringBuilder b = new StringBuilder( 512 );
        
        int space_cnt = 0;
        
        for( int i = 0; len == -1 || i < len; ++ i )
        {
            int c = istream.read();
            if( c == ' ' )
            {
                space_cnt++;
            } else if( c == '\r' )
            {
                --i;
                continue;
            } else if( c == '\n' )
            {   
                line_no ++;
                //get rid of trailing spaces;
                if( space_cnt > 0 )
                {
                    b.delete( b.length() - space_cnt, b.length() );
                    i -= space_cnt;
                    space_cnt -= space_cnt;
                }
                
                //return if no actual stop symbol and len is infinite
                if( -1 == stop && -1 == len ) 
                    break;
                --i;
                continue;
            } else if( c == -1 )
            {
                //TODO: should remove trailing spaces?
                if( b.length() > 0 )
                    return b.toString();
                
                throw new EOFException();        
            } else
            {
                space_cnt -= space_cnt;
            }
            
            b.append( (char)c );
            if( c == stop )
                break;
        } 
        return b.toString();
    }
}
