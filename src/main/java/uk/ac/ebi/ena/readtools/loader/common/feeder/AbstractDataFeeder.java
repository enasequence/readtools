package uk.ac.ebi.ena.readtools.loader.common.feeder;

import java.io.BufferedInputStream;
import java.io.EOFException;
import java.io.InputStream;
import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.List;

import uk.ac.ebi.ena.readtools.loader.common.eater.DataEater;
import uk.ac.ebi.ena.readtools.loader.common.eater.DataEaterException;



public abstract class 
AbstractDataFeeder<T> extends Thread implements DataFeeder<T>
{
    protected List<Method> feedables = new ArrayList<Method>();
    protected Method       checker;
    protected InputStream  istream;
    protected boolean      is_ok = true;
    protected DataEater<T, ?> dataEater;
    protected Throwable    stored_exception;
    private long           field_feed_count;
    static final int       YIELD_CYCLES = 362;//16384;
    
    
    protected 
    AbstractDataFeeder( InputStream istream, 
                        Class<T>    sample ) throws DataFeederException, SecurityException, NoSuchMethodException
    {
        this.istream = new BufferedInputStream( istream, 1024 * 1024 );
        Class<?> c = sample;
        while ( c != Object.class )
        {
            for( Field f : c.getDeclaredFields() )
            {
               FeedableData a = f.getAnnotation( FeedableData.class );
               if( null != a )
                   feedables.add( sample.getDeclaredMethod( a.method(), InputStream.class ) );
            }
            
            if( null == checker )
                for( Method m : c.getDeclaredMethods() )
                {
                    FeedableDataChecker a = m.getAnnotation( FeedableDataChecker.class );
                    if( null != a )
                    {
                        checker = m;
                        break;
                    }
                }
            c = c.getSuperclass();
        }
        
        
        if( feedables.size() == 0 )
            throw new DataFeederException( -1, "Sample structure does not contains public members annotated with @FeedableData" );
        
    }
    
    
    //Re-implement to instantiate feedable type of yours 
    protected abstract T
    newFeedable();
    
    public long
    getFieldFeedCount()
    {
        return field_feed_count;
    }
    
    //Re-implement if you need special type of feeding
    public T
    feed() throws DataFeederException 
    {
        T object = newFeedable();
        boolean record_started = false;
        try
        {
            try
            {   
                for( Method m : feedables )
                {
                    m.invoke( object, new Object[] { istream } );
                    field_feed_count++;
                    record_started = true;
                }
            } catch( InvocationTargetException ite )
            {
                Throwable cause = ite.getCause();
                if( cause instanceof EOFException )
                {
                	if( record_started )
                    {
                        if( null != checker )
                        {
                            checker.invoke( object, (Object [])null );
                        } else
                        {
                        throw new DataFeederException( field_feed_count, "EOF while reading fields" );
                        }
                    }
                	throw new DataFeederEOFException( field_feed_count );
                } else if( cause instanceof DataFeederException )
                {
                	throw (DataFeederException) cause;
                	
                } else
                {
                	throw ite;
                }
            } 
            
            if( null != checker )
                checker.invoke( object, (Object [])null );
            
            return object;        

        } catch( IllegalAccessException | IllegalArgumentException | InvocationTargetException e)
        {
            throw new DataFeederPanicException( e );
        }
    }
 
    
    public AbstractDataFeeder<T>
    setEater( DataEater<T, ?> eater )
    {
        this.dataEater = eater;
        return this;
    }
    
    //TODO: scheduler should be fair and based on different principle
    public final void
    run()
    {
        try
        {
            int yield = YIELD_CYCLES;
            for( ; ; )
            {
                synchronized( dataEater )
                {
                    for( yield = YIELD_CYCLES; yield > 0; --yield )
                        dataEater.eat( feed() );
                }
                
                if( !dataEater.isOk() )
                    throw new DataFeederPanicException();
                
                Thread.sleep( 1 );
            }
            
        } catch( DataEaterException e )
        {
            //e.printStackTrace();
            this.stored_exception = e;
            is_ok = false;
        } catch( DataFeederEOFException e )
        {
            //e.printStackTrace();
        } catch( DataFeederPanicException e )
        {
            //e.printStackTrace();
            is_ok = false;
            this.stored_exception = e;
        } catch( Throwable t )
        {
            //t.printStackTrace();
            this.stored_exception = t;
            is_ok = false;
        }
    }
    
    
    public boolean
    isOk()
    {
        return is_ok;
    }
    
    
    
    public Throwable
    getStoredException()
    {
        return stored_exception;
    }
}
