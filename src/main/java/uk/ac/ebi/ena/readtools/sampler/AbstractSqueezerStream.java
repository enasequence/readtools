package uk.ac.ebi.ena.readtools.sampler;

import java.io.IOException;
import java.io.OutputStream;
import java.util.BitSet;

public abstract class 
AbstractSqueezerStream
{
    OutputStream ostream;
    BitSet       remains;
    
    public void
    setOutputStream( OutputStream ostream )
    {
        this.ostream = ostream;
    }
    
    
    BitSet
    writeSet( BitSet set ) throws IOException 
    {
        BitSet remains = new BitSet( 0 );
        int mod = set.size() % 8;
        if( 0 != mod )
        {
            remains = set.get( set.size() - mod, set.size() );
        }
        
        int upto = set.size() - mod;
        int value = 0;
        for( int i = 0; i < upto; ++i )
        {
            if( set.get( i ) )
                value = value | 1 << ( i % 8 );
         
            if( 0 == i % 8 )
            {
                ostream.write( value );
                value ^=value;
            }
        }    
        return remains;
    }
    
    public abstract BitSet
    squeezeMethod( BitSet remains, byte bytes[] );

    
    public void
    squeeze( byte bytes[] ) throws IOException
    {
        BitSet b = squeezeMethod( remains, bytes );
        remains = writeSet( b );
    }
    
    
    public void
    finish()
    {
        
    }
}
