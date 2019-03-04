package uk.ac.ebi.ena.sampler.bits;

public interface 
ConstantLengthDataDictionary<T>
{
    public int getIncomingBlockLengthBytes();   //bytes!
    public int getTranslatedLength(); //bits!
    public int translate( int value );
    public int translateBack( int value );
    public int getDictionarySize();
    
}
