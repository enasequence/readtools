package uk.ac.ebi.ena.frankenstein.loader.common.eater;

public class 
DataEaterException extends Exception
{
    private static final long serialVersionUID = 1L;

    public 
    DataEaterException( String value )
    {
        super( value );
    }

    public 
    DataEaterException()
    {
        super();
    }

    
    public 
    DataEaterException( Throwable cause )
    {
        super( cause );
    }

}
