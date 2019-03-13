package uk.ac.ebi.ena.readtools.loader.fastq;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.SynchronousQueue;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicReference;

import uk.ac.ebi.ena.readtools.loader.common.FileCompression;
import uk.ac.ebi.ena.readtools.loader.common.QualityNormalizer;
import uk.ac.ebi.ena.readtools.loader.common.eater.DataEater;
import uk.ac.ebi.ena.readtools.loader.common.eater.DataEaterException;
import uk.ac.ebi.ena.readtools.loader.common.feeder.AbstractDataFeeder;
import uk.ac.ebi.ena.readtools.loader.common.feeder.DataFeederException;
import uk.ac.ebi.ena.readtools.loader.fastq.DataSpot.DataSpotParams;
import uk.ac.ebi.ena.readtools.loader.fastq.IlluminaIterativeEater.READ_TYPE;

public class 
IlluminaIterativeEaterIterator implements Iterator<IlluminaSpot>, DataEater<IlluminaSpot, Object>
{
    private static final long CYCLE_TIMEFRAME = 0L;
    private BlockingQueue<IlluminaSpot>   queue = new SynchronousQueue<IlluminaSpot>();
    private AtomicReference<IlluminaSpot> current_element = new AtomicReference<IlluminaSpot>();
    private AtomicBoolean was_cascade_errors = new AtomicBoolean( false );
    private Exception storedException;    
    

    IlluminaIterativeEaterIterator( File tmp_folder, 
                            int spill_page_size, //only for paired 
                            READ_TYPE read_type, 
                            File[] files, 
                            final QualityNormalizer normalizers[]  ) throws SecurityException, FileNotFoundException, DataFeederException, NoSuchMethodException, IOException
    {
        DataEater<DataSpot, IlluminaSpot> eater = null;
        
        
        switch( read_type )
        {
        
        case SINGLE:
            eater = (DataEater<DataSpot, IlluminaSpot>)new IlluminaSingleDataEater();
            break;
            
        case PAIRED:
            eater = (DataEater<DataSpot, IlluminaSpot>) new IlluminaPairedDataEater( tmp_folder, spill_page_size );
            break;
            
        default:
            throw new UnsupportedOperationException();

        }
        
        //IterativeEater vdb2  = new IterativeEater();//type, read_type, p.data_folder_name, p.use_md5 );
        
        eater.setEater( this );
        
        ArrayList<AbstractDataFeeder<?>> feeders = new ArrayList<AbstractDataFeeder<?>>();
        

        int attr = 1;
        
        for( File file: files )
        {   
            final String default_attr = Integer.toString( attr ++ );
            final int nindex = normalizers.length == files.length ? attr - 2 : 0;
            
            AbstractDataFeeder<?> feeder = 
            new AbstractDataFeeder<DataSpot>( FileCompression.open( file ), DataSpot.class ) 
            {
                final DataSpotParams params = DataSpot.defaultParams();
                
                @Override
                protected DataSpot 
                newFeedable()
                {
                    return new DataSpot( normalizers[ nindex ], default_attr, params );
                }
            }.setEater( eater );
            
            feeder.setName( file.getPath() );
            feeders.add( feeder );
            feeder.start();
        }
        
        new Thread( lifecycle( feeders, eater ), "lifecycle" ).start();
    }
    
    
    private Runnable
    lifecycle( final ArrayList<AbstractDataFeeder<?>> feeders, 
               final DataEater<?, ?> eater_root )
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
                                for( AbstractDataFeeder<?> feeder : feeders )
                                {
                                    if( feeder.isAlive() )
                                    {
                                        try
                                        {
                                            feeder.join( CYCLE_TIMEFRAME );
                                        } catch( InterruptedException ie )
                                        {
                                            again = true;
                                            System.out.printf( "%s was interrupted\n", feeder.getName() );
                                        }
                                    }else if( !feeder.isOk() )
                                    {
                                        throw new DataFeederException( feeder.getStoredException() );
                                    }
                                }
                            }while( again );
                           
                            for( AbstractDataFeeder<?> feeder : feeders )
                            {
                                if( !feeder.isOk() )
                                {
                                    throw new DataFeederException( feeder.getStoredException() );
                                }
                            }
                    
                            eater_root.cascadeErrors();

                        }catch( Exception dfe )
                        {
                            storedException = dfe;
                        }
                    }
               };
    }
    
    
    @Override
    public void 
    cascadeErrors() throws DataEaterException
    {
        was_cascade_errors.lazySet( true );
    }

    
    @Override
    public void 
    eat( IlluminaSpot object ) throws DataEaterException
    {
        try
        {
            queue.put( object );
            
        } catch( InterruptedException e )
        {
            e.printStackTrace();
        }
    }

    
    @Override
    public void 
    setEater( DataEater<Object, ?> dataEater )
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
    public IlluminaSpot 
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
