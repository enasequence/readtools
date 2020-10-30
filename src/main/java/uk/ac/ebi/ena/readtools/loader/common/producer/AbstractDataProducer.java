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

import uk.ac.ebi.ena.readtools.loader.common.consumer.Spot;
import uk.ac.ebi.ena.readtools.loader.common.consumer.DataConsumer;
import uk.ac.ebi.ena.readtools.loader.common.consumer.DataConsumerException;

import java.io.BufferedInputStream;
import java.io.EOFException;
import java.io.InputStream;

public abstract class
AbstractDataProducer<T extends Spot> extends Thread implements DataProducer<T> {
    private static final int YIELD_CYCLES = 362;//16384;

    private volatile long recordCount = 0, totalBaseCount = 0;

    protected final InputStream  istream;

    protected volatile DataConsumer<T, ?> dataConsumer;
    protected volatile boolean      is_ok = true;
    protected volatile Throwable    stored_exception;
    
    protected AbstractDataProducer(InputStream istream) {
        this.istream = new BufferedInputStream( istream, 1024 * 1024 );
    }

    /**
     * Get the total number of records that were read.
     *
     * @return
     */
    public long getRecordCount() {
        return recordCount;
    }

    /**
     * Get the total number of bases that were read.
     *
     * @return
     */
    public long getTotalBaseCount() {
        return totalBaseCount;
    }
    
    public void setConsumer(DataConsumer<T, ?> consumer) {
        this.dataConsumer = consumer;
    }
    
    public boolean isOk() {
        return is_ok;
    }
    
    public Throwable getStoredException() {
        return stored_exception;
    }

    //TODO: scheduler should be fair and based on different principle
    public final void run() {
        try {
            begin();

            for( ; ; ) {
                synchronized(dataConsumer) {
                    for(int yield = YIELD_CYCLES; yield > 0; --yield )
                        dataConsumer.consume( produce() );
                }

                if( !dataConsumer.isOk() )
                    throw new DataProducerPanicException();

                Thread.sleep( 1 );
            }
        } catch( DataConsumerException e ) {
            //e.printStackTrace();
            this.stored_exception = e;
            is_ok = false;
        } catch( DataProducerEOFException e ) {
            //e.printStackTrace();
        } catch( DataProducerPanicException e ) {
            //e.printStackTrace();
            is_ok = false;
            this.stored_exception = e;
        } catch( Throwable t ) {
            //t.printStackTrace();
            this.stored_exception = t;
            is_ok = false;
        } finally {
            end();
        }
    }

    protected void begin() {}

    protected void end() {}

    //Re-implement if you need special type of feeding
    private T produce() {
        T spot = null;

        try {
            spot = produce(istream);
            ++recordCount;
            totalBaseCount += spot.getBaseCount();


            return spot;
        } catch( EOFException e ){
            throw new DataProducerEOFException(recordCount);
        } catch( DataProducerException e ){
            throw e;
        } catch( Throwable cause ) {
            throw new DataProducerException(cause);
        }
    }
}
