package uk.ac.ebi.ena.sampler;

public class 
BaseSqueezer
{
    BaseSqueezer()
    {
        
    }
    
    
    byte[]
    squeeze( byte value[] )
    {
        byte[] result = new byte[ value.length / 2 + value.length % 2 ];
        int r  = 0;
        for( int i = 0; i < value.length; i += 2 )
            result[ r++ ] = (byte) squeezePair( value[ i ], 
                                                i + 1 < value.length ? value[ i + 1 ] : (byte)'!' );
        
        return result;
    }
    
    
    int
    squeezePair( byte one, byte two )
    {
       return ( translate( one ) << 4 ) | ( 0xF & translate( two ) ); 
    }
    

    int
    translate( int value )
    {
        switch( value )
        {
        case 'A': return 0x00;
        case 'C': return 0x01;
        case 'G': return 0x02;
        case 'T': return 0x03;
        case 'N': return 0x04;
        
        default: 
            throw new RuntimeException( String.format( "byte %s could not be translated", value ) );
        }

    }
    
    
}
