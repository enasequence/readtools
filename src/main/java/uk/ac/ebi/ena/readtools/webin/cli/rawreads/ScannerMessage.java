package uk.ac.ebi.ena.readtools.webin.cli.rawreads;

public interface 
ScannerMessage
{
	String    getMessage();
	String    getOrigin();
	Throwable getThrowable();
	
	abstract static class 
	AbstractMessage implements ScannerMessage
	{
		@Override public String 
		toString()
		{
			return String.format( "%s: %s%s%s", getClass().getSimpleName(), 
					              getMessage(), 
					              null == getOrigin() ? "" : "<" + getOrigin() + ">", 
					              null == getThrowable() ? "" : getThrowable() );
		}
	}
	
	
	public static class 
	ScannerInfoMessage extends AbstractMessage implements ScannerMessage 
	{
		private Throwable throwable;
		private String    message;
		private String    origin;
		
		
		public
		ScannerInfoMessage( String m )
		{
			this( m, null );
		}
	
		
		public
		ScannerInfoMessage( String m, String o )
		{
			this.message = m;
			this.origin  = o;
		}
		
		
		@Override public String
		getMessage()
		{
			return message;
		}
		
	
		@Override public String
		getOrigin()
		{
			return origin;
		}
	
		
		@Override public Throwable
		getThrowable()
		{
			return this.throwable;
		}
	}
	
	
	public static class 
	ScannerErrorMessage extends AbstractMessage implements ScannerMessage
	{
		Throwable throwable;
		String    message;
		String    origin;
		
		
		public
		ScannerErrorMessage( String m )
		{
			this( null, m, null );
		}
		
		
		public
		ScannerErrorMessage( Throwable t, String m, String o )
		{
			this.throwable = t;
			this.message   = m;
			this.origin    = o;
		}
		
		
		@Override public String
		getMessage()
		{
			return message;
		}
		

		@Override public String
		getOrigin()
		{
			return origin;
		}

		
		@Override public Throwable
		getThrowable()
		{
			return this.throwable;
		}
	}

}