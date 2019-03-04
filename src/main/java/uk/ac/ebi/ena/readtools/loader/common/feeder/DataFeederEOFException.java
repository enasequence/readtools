package uk.ac.ebi.ena.readtools.loader.common.feeder;

public class 
DataFeederEOFException extends DataFeederException
{
    private static final long serialVersionUID = 1L;
    private final String thread_name;
    
    
    public 
    DataFeederEOFException( String value )
    {
        super( -1, value );
        thread_name = Thread.currentThread().getName();
    }

    public 
    DataFeederEOFException()
    {
        super( -1 );
        thread_name = Thread.currentThread().getName();
    }
    
    public 
    DataFeederEOFException( Throwable cause )
    {
        super( cause );
        thread_name = Thread.currentThread().getName();
    }
    
    
    public String
    toString()
    {
        return String.format( "%s: %s", thread_name, super.toString() ); 
    }
    
}
