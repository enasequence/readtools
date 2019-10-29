package uk.ac.ebi.ena.readtools.loader.common.feeder;

public class 
DataFeederPanicException extends DataFeederException
{
    private static final long serialVersionUID = 1L;
    private final String thread_name;
    
    
    public 
    DataFeederPanicException( String value )
    {
        super( -1, value );
        thread_name = Thread.currentThread().getName();
    }

    
    public 
    DataFeederPanicException()
    {
        super( -1 );
        thread_name = Thread.currentThread().getName();
    }
    
    
    public 
    DataFeederPanicException( Throwable cause )
    {
        super( cause );
        thread_name = Thread.currentThread().getName();
    }
    
    
    public 
    DataFeederPanicException( String message, Throwable cause )
    {
        super( message, cause );
        thread_name = Thread.currentThread().getName();
    }
    
    
    public String
    toString()
    {
        return String.format( "%s: %s", thread_name, super.toString() ); 
    }
    
}
