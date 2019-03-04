package uk.ac.ebi.ena.frankenstein.loader.cg.reads;

import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;

import uk.ac.ebi.ena.frankenstein.loader.common.FileCompression;
import uk.ac.ebi.ena.frankenstein.loader.common.eater.DataEater;
import uk.ac.ebi.ena.frankenstein.loader.common.eater.DataEaterException;
import uk.ac.ebi.ena.frankenstein.loader.common.eater.PrintDataEater;
import uk.ac.ebi.ena.frankenstein.loader.common.feeder.AbstractDataFeeder;
import uk.ac.ebi.ena.frankenstein.loader.common.feeder.DataFeederException;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;

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
          DataEater<CGSpot, T> next_eater ) throws DataFeederException, IOException, SecurityException, NoSuchMethodException, DataEaterException
    {
        ArrayList<AbstractDataFeeder<?>> feeders = new ArrayList<AbstractDataFeeder<?>>();
        
        final ReadHeader read_header = ReadHeader( read_stream );
        final ReadHeader map_header  = ReadHeader( map_stream );
        if( !read_header.getFileProrerty( ReadHeader.PropertyName.BATCH_FILE_NUMBER ).equals( map_header.getFileProrerty( ReadHeader.PropertyName.BATCH_FILE_NUMBER ) ) )
                throw new IOException( "Batches did not match!" );
        
        AbstractDataFeeder<CompleteGenomicsMap3> mf = new AbstractDataFeeder<CompleteGenomicsMap3>( map_stream, CompleteGenomicsMap3.class ) 
        {

            @Override
            protected CompleteGenomicsMap3 
            newFeedable()
            {
                return new CompleteGenomicsMap3( map_header );
            }
        };

        AbstractDataFeeder<CompleteGenomicsRead> rf = new AbstractDataFeeder<CompleteGenomicsRead>( read_stream, CompleteGenomicsRead.class ) 
        {

            @Override
            protected CompleteGenomicsRead 
            newFeedable()
            {
                return new CompleteGenomicsRead( read_header );
            }
        };
        

        
        DataEater<? extends CompleteGenomicsBase, CGSpot> eater = new CompleteGenomicsPairedEater();
        eater.setEater( next_eater );        
        mf.setEater( (DataEater<CompleteGenomicsMap3, ?>) eater );
        feeders.add( mf );
        mf.start();
        
        rf.setEater( (DataEater<CompleteGenomicsRead, ?>) eater );
        feeders.add( rf );
        rf.start();

        boolean again = false;
        do
        {
            for( AbstractDataFeeder<?> feeder : feeders )
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
    
            read( read_stream, map_stream, (DataEater<CGSpot, ?>)new PrintDataEater() );
            
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