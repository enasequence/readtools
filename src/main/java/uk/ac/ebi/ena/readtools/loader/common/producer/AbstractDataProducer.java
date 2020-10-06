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

import uk.ac.ebi.ena.readtools.loader.common.consumer.DataConsumer;
import uk.ac.ebi.ena.readtools.loader.common.consumer.DataConsumerException;



public abstract class
AbstractDataProducer<T extends DataProducible> extends Thread implements DataProducer<T>
{
    protected InputStream  istream;
    protected boolean      is_ok = true;
    protected DataConsumer<T, ?> dataConsumer;
    protected Throwable    stored_exception;
    private long readRecordCount;
    static final int       YIELD_CYCLES = 362;//16384;
    
    
    protected AbstractDataProducer(InputStream istream)
    {
        this.istream = new BufferedInputStream( istream, 1024 * 1024 );
    }
    
    
    //Re-implement to instantiate feedable type of yours 
    protected abstract T
    newProducible();

    /**
     * Get the total number of records that were read.
     *
     * @return
     */
    public long
    getReadRecordCount()
    {
        return readRecordCount;
    }

    //Re-implement if you need special type of feeding
    public T
    produce() throws DataProducerException
    {
        T object = newProducible();

        try
        {
            try
            {
                object.read(istream);
                ++readRecordCount;
            } catch( EOFException e ){
                throw new DataProducerEOFException(readRecordCount);
            } catch( DataProducerException e ){
                throw e;
            } catch( Throwable cause ) {
                throw new DataProducerException(cause);
            }

            return object;

        } catch( IllegalArgumentException e) {
            throw new DataProducerPanicException( String.valueOf( object ), e );
        }
    }
 
    
    public void setConsumer(DataConsumer<T, ?> consumer)
    {
        this.dataConsumer = consumer;
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