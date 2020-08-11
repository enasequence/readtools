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
package uk.ac.ebi.ena.readtools.sampler.bits;

import java.io.IOException;
import java.io.OutputStream;


public class
DictionaryOutputStream extends BitOutputStream
{
    final static int BYTE_LENGTH = 8; 
    ConstantLengthDataDictionary<?> dict;
    int block_size = -1;
    
    
    public 
    DictionaryOutputStream( OutputStream   out, 
                            ConstantLengthDataDictionary<?> dict ) throws IOException
    {
        super( out );
        this.dict = dict;
        block_size = dict.getIncomingBlockLengthBytes();
        if( block_size > 1 )
            throw new IOException( "Datablocks more than 1 byte are not supported" );
    }

    public void 
    write( byte[] bytes ) throws IOException
    {
        if( ( bytes.length % block_size ) != 0 )
            throw new IOException( "Incoming data size not fits into block size" );
        
        for( byte b : bytes )
            super.write( dict.translate( b ), dict.getTranslatedLength() );
    }
    
    
}
