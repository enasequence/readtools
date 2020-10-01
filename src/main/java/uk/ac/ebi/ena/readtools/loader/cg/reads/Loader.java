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
package uk.ac.ebi.ena.readtools.loader.cg.reads;

import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;

import uk.ac.ebi.ena.readtools.loader.common.FileCompression;
import uk.ac.ebi.ena.readtools.loader.common.eater.DataConsumer;
import uk.ac.ebi.ena.readtools.loader.common.eater.DataConsumerException;
import uk.ac.ebi.ena.readtools.loader.common.eater.PrintDataEater;
import uk.ac.ebi.ena.readtools.loader.common.feeder.AbstractDataProducer;
import uk.ac.ebi.ena.readtools.loader.common.feeder.DataProducerException;

public class 
Loader
{
    static ReadHeader
    ReadHeader( InputStream is ) throws IOException
    {
        ReadHeader header = new ReadHeader( is );
        header.read();
        System.out.println( header.toString() );
        return header;
    }
    
    
    public static <T> void 
    read( InputStream read_stream, 
          InputStream map_stream, 
          DataConsumer<CGSpot, T> next_eater ) throws DataProducerException, IOException, SecurityException, NoSuchMethodException, DataConsumerException
    {
        ArrayList<AbstractDataProducer<?>> feeders = new ArrayList<AbstractDataProducer<?>>();
        
        final ReadHeader read_header = ReadHeader( read_stream );
        final ReadHeader map_header  = ReadHeader( map_stream );
        if( !read_header.getFileProrerty( ReadHeader.PropertyName.BATCH_FILE_NUMBER ).equals( map_header.getFileProrerty( ReadHeader.PropertyName.BATCH_FILE_NUMBER ) ) )
                throw new IOException( "Batches did not match!" );
        
        AbstractDataProducer<CompleteGenomicsMap3> mf = new AbstractDataProducer<CompleteGenomicsMap3>( map_stream, CompleteGenomicsMap3.class )
        {

            @Override
            protected CompleteGenomicsMap3
            newProducible()
            {
                return new CompleteGenomicsMap3( map_header );
            }
        };

        AbstractDataProducer<CompleteGenomicsRead> rf = new AbstractDataProducer<CompleteGenomicsRead>( read_stream, CompleteGenomicsRead.class )
        {

            @Override
            protected CompleteGenomicsRead
            newProducible()
            {
                return new CompleteGenomicsRead( read_header );
            }
        };
        

        
        DataConsumer<? extends CompleteGenomicsBase, CGSpot> eater = new CompleteGenomicsPairedEater();
        eater.setConsumer( next_eater );
        mf.setConsumer( (DataConsumer<CompleteGenomicsMap3, ?>) eater );
        feeders.add( mf );
        mf.start();
        
        rf.setConsumer( (DataConsumer<CompleteGenomicsRead, ?>) eater );
        feeders.add( rf );
        rf.start();

        boolean again = false;
        do
        {
            for( AbstractDataProducer<?> feeder : feeders )
            {
                if( feeder.isAlive() )
                {
                    try
                    {
                        feeder.join();
                    } catch( InterruptedException ie )
                    {
                        again = true;
                    }
                }else if( !feeder.isOk() )
                {
                    throw new DataProducerException( feeder.getStoredException() );
                }
            }
        }while( again );
       
        for( AbstractDataProducer<?> feeder : feeders )
        {
            if( !feeder.isOk() )
            {
                throw new DataProducerException( feeder.getStoredException() );
            }
        }
    
    
        eater.cascadeErrors();
    }


    public static void
    main( String args[] ) throws Exception
    {
        Params params = new Params();
        JCommander jc = new JCommander( params );
        try
        {
            jc.parse( args );
            
        }catch( Exception e )
        {
            jc.usage();
            System.exit( 1 );
            
        }
    
        InputStream read_stream = null;
        InputStream map_stream  = null;
        try
        {
            read_stream = params.bz2 ? FileCompression.BZ2.open( params.read, false ) : FileCompression.NONE.open( params.read, false );
            map_stream  = params.bz2 ? FileCompression.BZ2.open( params.map, false ) : FileCompression.NONE.open( params.map, false );
    
            read( read_stream, map_stream, (DataConsumer<CGSpot, ?>)new PrintDataEater() );
            
        } finally
        {
            if( null != read_stream )
                read_stream.close();
    
            if( null != map_stream )
                map_stream.close();
        }
    }
    
    
    static class 
    Params
    {
        @Parameter( names = "--read", required = true )
        String read;
        
        @Parameter( names = "--map", required = true )
        String map;
        
        @Parameter( names = "--bz2" )
        boolean bz2;
        
    }

}