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
package uk.ac.ebi.ena.readtools.sampler.bits;

import java.io.IOException;


public class
DictionaryInputStream extends BitInputStream
{
    ConstantLengthDataDictionary<?> dict;
    int block_size  = -1;
    int tblock_size = -1;
    
    public 
    DictionaryInputStream( BitInputStream in, 
                           ConstantLengthDataDictionary<?> dict ) throws IOException
    {
        super( in );
        this.dict = dict;
        block_size = dict.getIncomingBlockLengthBytes();
        if( block_size > 1 )
            throw new IOException( "Datablocks more than 1 byte are not supported" );
        
        tblock_size = dict.getTranslatedLength();
    }
    
    public int 
    read() throws IOException
    {
        return dict.translateBack( read( tblock_size ) );
    }
    
    
}
