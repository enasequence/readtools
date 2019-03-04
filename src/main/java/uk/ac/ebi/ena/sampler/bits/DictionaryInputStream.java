package uk.ac.ebi.ena.sampler.bits;

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
