/*
* Copyright 2010-2020 EMBL - European Bioinformatics Institute
* Licensed under the Apache License, Version 2.0 (the "License"); you may not use this
* file except in compliance with the License. You may obtain a copy of the License at
* http://www.apache.org/licenses/LICENSE-2.0
* Unless required by applicable law or agreed to in writing, software distributed under the
* License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
* CONDITIONS OF ANY KIND, either express or implied. See the License for the
* specific language governing permissions and limitations under the License.
*/
package uk.ac.ebi.ena.readtools.loader.common.producer;

import java.io.BufferedInputStream;
import java.io.EOFException;
import java.io.InputStream;
import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.List;

import uk.ac.ebi.ena.readtools.loader.common.consumer.DataConsumable;
import uk.ac.ebi.ena.readtools.loader.common.consumer.DataConsumer;
import uk.ac.ebi.ena.readtools.loader.common.consumer.DataConsumerException;



public abstract class
AbstractDataProducer<T extends DataConsumable> extends Thread implements DataProducer<T>
{
    protected List<Method> producibles = new ArrayList<Method>();
    protected Method       checker;
    protected InputStream  istream;
    protected boolean      is_ok = true;
    protected DataConsumer<T, ?> dataConsumer;
    protected Throwable    stored_exception;
    private long           field_feed_count;
    static final int       YIELD_CYCLES = 362;//16384;
    
    
    protected AbstractDataProducer(InputStream istream,
                                   Class<T>    sample ) throws DataProducerException, SecurityException, NoSuchMethodException
    {
        this.istream = new BufferedInputStream( istream, 1024 * 1024 );
        Class<?> c = sample;
        while ( c != Object.class )
        {
            for( Field f : c.getDeclaredFields() )
            {
               ProducibleData a = f.getAnnotation( ProducibleData.class );
               if( null != a )
                   producibles.add( sample.getDeclaredMethod( a.method(), InputStream.class ) );
            }
            
            if( null == checker )
                for( Method m : c.getDeclaredMethods() )
                {
                    ProducibleDataChecker a = m.getAnnotation( ProducibleDataChecker.class );
                    if( null != a )
                    {
                        checker = m;
                        break;
                    }
                }
            c = c.getSuperclass();
        }
        
        
        if( producibles.size() == 0 )
            throw new DataProducerException( -1, "Sample structure does not contains public members annotated with @FeedableData" );
        
    }
    
    
    //Re-implement to instantiate feedable type of yours 
    protected abstract T
    newProducible();
    
    public long
    getFieldFeedCount()
    {
        return field_feed_count;
    }
    
    //Re-implement if you need special type of feeding
    public T
    produce() throws DataProducerException
    {
        T object = newProducible();
        boolean record_started = false;
        try
        {
            try
            {   
                for( Method m : producibles)
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
                            throw new DataProducerException( field_feed_count, "EOF while reading fields" );
                        }
                    }
                	throw new DataProducerEOFException( field_feed_count );
                } else if( cause instanceof DataProducerException)
                {
                	throw (DataProducerException) cause;
                	
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
            throw new DataProducerPanicException( String.valueOf( object ), e );
        }
    }
 
    
    public AbstractDataProducer<T>
    setConsumer(DataConsumer<T, ?> consumer)
    {
        this.dataConsumer = consumer;
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
                synchronized(dataConsumer)
                {
                    for( yield = YIELD_CYCLES; yield > 0; --yield )
                        dataConsumer.consume( produce() );
                }
                
                if( !dataConsumer.isOk() )
                    throw new DataProducerPanicException();
                
                Thread.sleep( 1 );
            }
            
        } catch( DataConsumerException e )
        {
            //e.printStackTrace();
            this.stored_exception = e;
            is_ok = false;
        } catch( DataProducerEOFException e )
        {
            //e.printStackTrace();
        } catch( DataProducerPanicException e )
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
