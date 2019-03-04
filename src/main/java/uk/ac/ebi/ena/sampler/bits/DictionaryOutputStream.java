package uk.ac.ebi.ena.sampler.bits;

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
