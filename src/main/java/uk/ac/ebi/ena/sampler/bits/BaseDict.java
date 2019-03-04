package uk.ac.ebi.ena.sampler.bits;


public class 
BaseDict implements ConstantLengthDataDictionary<Character>
{

    @Override
    public int 
    getIncomingBlockLengthBytes()
    {
        return 1;
    }

    @Override
    public int 
    getTranslatedLength()
    {
        return 3;
    }
/*
    7   111
    a   4   010
    e   4   000
    f   3   1101
    h   2   1010
    
    i   2   1000
    m   2   0111
    n   2   0010
    s   2   1011
    t   2   0110
*/    
    @Override
    public int 
    translate( int value )
    {
        switch( value )
        {
        case 'A': return 0x8;
        case 'T': return 0x7;
        case 'G': return 0x2;
        case 'C': return 0x11;
        case 'N': return 0x6;
        
        default:
            throw new RuntimeException();
        }
    }


    @Override
    public int
    translateBack( int value )
    {
        switch( value )
        {
        case 0x0: return 'A';
        case 0x1: return 'C';
        case 0x2: return 'G';
        case 0x3: return 'T';
        case 0x4: return 'N';
        default:
            throw new RuntimeException();
        }
    }
 
    
    
    @Override
    public int 
    getDictionarySize()
    {
    
        return 5;
    }
    
}
