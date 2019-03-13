package uk.ac.ebi.ena.readtools.loader.common.feeder;

public class 
DataFeederException extends Exception
{
    private static final long serialVersionUID = 1L;
    private long   line_no = -1;
    private final  String thread_name = Thread.currentThread().getName();
    
    public 
    DataFeederException( long line_no, String value )
    {
        super( value );
        this.line_no = line_no;
    }

    public 
    DataFeederException( long line_no )
    {
        super();
        this.line_no = line_no;
    }

    
    public 
    DataFeederException( long line_no, Throwable cause )
    {
        super( cause );
        this.line_no = line_no;
    }


    public 
    DataFeederException( Throwable cause )
    {
        super( cause );
    }
   
    
    public String
    toString()
    {
        return line_no < 0 ? super.toString() 
                           : String.format( "%s:%d %s", thread_name, line_no, super.toString() );
    }
    
    
    public long
    getLineNo()
    {
    	return this.line_no;
    }
}
