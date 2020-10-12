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
package uk.ac.ebi.ena.readtools.loader.fastq;

import uk.ac.ebi.ena.readtools.common.reads.QualityNormalizer;
import uk.ac.ebi.ena.readtools.loader.common.FileCompression;
import uk.ac.ebi.ena.readtools.loader.common.consumer.Spot;
import uk.ac.ebi.ena.readtools.loader.common.consumer.DataConsumer;
import uk.ac.ebi.ena.readtools.loader.common.consumer.DataConsumerException;
import uk.ac.ebi.ena.readtools.loader.common.producer.AbstractDataProducer;
import uk.ac.ebi.ena.readtools.loader.common.producer.DataProducerException;
import uk.ac.ebi.ena.readtools.loader.common.producer.DataSpotProducer;
import uk.ac.ebi.ena.readtools.loader.fastq.FastqIterativeConsumer.READ_TYPE;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.SynchronousQueue;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicReference;

public class
FastqIterativeConsumerIterator implements Iterator<FastqSpot>, DataConsumer<FastqSpot, Spot>
{
    private static final long CYCLE_TIMEFRAME = 0L;
    private BlockingQueue<FastqSpot>   queue = new SynchronousQueue<FastqSpot>();
    private AtomicReference<FastqSpot> current_element = new AtomicReference<FastqSpot>();
    private AtomicBoolean was_cascade_errors = new AtomicBoolean( false );
    private Exception storedException;    
    

    FastqIterativeConsumerIterator(File tmp_folder,
                                   int spill_page_size, //only for paired
                                   READ_TYPE read_type,
                                   File[] files,
                                   final QualityNormalizer normalizers[]  ) throws SecurityException, IOException
    {
        DataConsumer<DataSpot, FastqSpot> consumer = null;
        
        
        switch( read_type )
        {
        
        case SINGLE:
            consumer = new SingleFastqConsumer();
            break;
            
        case PAIRED:
            consumer = new PairedFastqConsumer( tmp_folder, spill_page_size );
            break;
            
        default:
            throw new UnsupportedOperationException();

        }

        consumer.setConsumer( this );
        
        ArrayList<AbstractDataProducer<?>> producers = new ArrayList<>();
        

        int attr = 1;
        
        for( File file: files )
        {   
            final String default_attr = Integer.toString( attr ++ );
            final int nindex = normalizers.length == files.length ? attr - 2 : 0;

            DataSpotProducer producer = new DataSpotProducer( FileCompression.open( file ), normalizers[ nindex ], default_attr );
            producer.setConsumer( consumer );
            producer.setName( file.getPath() );
            producers.add( producer );
            producer.start();
        }
        
        new Thread( lifecycle( producers, consumer ), "lifecycle" ).start();
    }
    
    
    private Runnable
    lifecycle( final ArrayList<AbstractDataProducer<?>> producers,
               final DataConsumer<?, ?> consumer_root )
    {
        return new Runnable()
               {
                    public void 
                    run()
                    {
                        try
                        {
                            boolean again = false;
                            do
                            {
                                for( AbstractDataProducer<?> producer : producers )
                                {
                                    if( producer.isAlive() )
                                    {
                                        try
                                        {
                                            producer.join( CYCLE_TIMEFRAME );
                                        } catch( InterruptedException ie )
                                        {
                                            again = true;
                                            System.out.printf( "%s was interrupted\n", producer.getName() );
                                        }
                                    }else if( !producer.isOk() )
                                    {
                                        throw new DataProducerException( producer.getStoredException() );
                                    }
                                }
                            }while( again );
                           
                            for( AbstractDataProducer<?> producer : producers )
                            {
                                if( !producer.isOk() )
                                {
                                    throw new DataProducerException( producer.getStoredException() );
                                }
                            }
                    
                            consumer_root.cascadeErrors();

                        }catch( Exception dfe )
                        {
                            storedException = dfe;
                        }
                    }
               };
    }
    
    
    @Override
    public void 
    cascadeErrors() throws DataConsumerException
    {
        was_cascade_errors.lazySet( true );
    }

    
    @Override
    public void
    consume(FastqSpot spot) throws DataConsumerException
    {
        try
        {
            queue.put(spot);
            
        } catch( InterruptedException e )
        {
            e.printStackTrace();
        }
    }

    
    @Override
    public void
    setConsumer(DataConsumer<Spot, ?> dataConsumer)
    {
        throw new RuntimeException( "N07 iMPl3m3nt3D" );
    }

    
    @Override
    public boolean 
    hasNext()
    {
        if( null != storedException )
            throw new RuntimeException( storedException );
        
        try
        {
            while( !was_cascade_errors.get() && false == current_element.compareAndSet( null, queue.poll( 1L, TimeUnit.SECONDS ) ) );
        } catch( InterruptedException e )
        {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
        return null != current_element.get();
    }
    
    
    
    @Override
    public FastqSpot
    next()
    {
        if( null != storedException )
            throw new RuntimeException( storedException );

        return current_element.getAndSet( null );
    }

    
    @Override
    public void 
    remove()
    {
        throw new RuntimeException( "N07 iMPl3m3nt3D" );
    }

    
    public Throwable
    getStoredException()
    {
        return storedException;
    }


    @Override
    public boolean 
    isOk()
    {
        return true;
    }
}
