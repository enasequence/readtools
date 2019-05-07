package uk.ac.ebi.ena.readtools.loader.cg.reads;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;

import org.junit.After;
import org.junit.Before;

import uk.ac.ebi.ena.readtools.loader.common.FileCompression;
import uk.ac.ebi.ena.readtools.loader.common.eater.DataEater;
import uk.ac.ebi.ena.readtools.loader.common.eater.DataEaterException;
import uk.ac.ebi.ena.readtools.loader.common.eater.PrintDataEater;
import uk.ac.ebi.ena.readtools.loader.common.feeder.AbstractDataFeeder;
import uk.ac.ebi.ena.readtools.loader.common.feeder.DataFeederException;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;

public class 
TestEater
{
    @Before
    public void
    init()
    {
       
    }
    
    @After
    public void
    unwind()
    {
        
    }
    

    boolean 
    read( InputStream read_stream, 
          InputStream map_stream ) throws SecurityException, DataFeederException, NoSuchMethodException, IOException, DataEaterException
    {
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
        

        
        DataEater<? extends CompleteGenomicsBase, ?> eater = new CompleteGenomicsPairedEater();
        eater.setEater( new PrintDataEater() );        
        mf.setEater( (DataEater<CompleteGenomicsMap3, ?>) eater );
        rf.setEater( (DataEater<CompleteGenomicsRead, ?>) eater );
        
        rf.run();
        mf.run();
        if( rf.isOk() && mf.isOk() )
            eater.cascadeErrors();
        else 
            return false;
        
        return true;
    }
    

    ReadHeader
    ReadHeader( InputStream is ) throws IOException
    {
        ReadHeader header = new ReadHeader( is );
        header.read();
        System.out.println( header.toString() );
        return header;
    }
    
    
    boolean 
    read( String read_resource, String map_resource ) throws Exception
    {
        InputStream read_stream = TestEater.class.getClassLoader().getResourceAsStream( read_resource );
        InputStream map_stream  = TestEater.class.getClassLoader().getResourceAsStream( map_resource );
        try
        {
            return read( read_stream, map_stream );
        } finally
        {
            read_stream.close();
            map_stream.close();
        }
    }
    
    
    boolean 
    read( File read_file, File map_file ) throws Exception
    {
        InputStream read_stream = new BufferedInputStream( new FileInputStream( read_file )  );
        InputStream map_stream  = new BufferedInputStream( new FileInputStream( map_file ) );
        try
        {
            return read( read_stream, map_stream );
        } finally
        {
            read_stream.close();
            map_stream.close();
        }
    }

    
    
    @org.junit.Test
    public void
    testCorrect() throws Exception
    {
     /*   
        if( !read( new File( "h:/mp3_schw3.fq" ), QualityNormalizer.ILLUMINA_1_3 ) )
            throw new Exception( "fail!" );
      */  
        if( !read( "cg_read.txt", "cg_map.txt" ) )
            throw new Exception( "fail!" );
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
        TestEater te = new TestEater();

        InputStream read_stream = null;
        InputStream map_stream  = null;
        try
        {
            read_stream = params.bz2 ? FileCompression.BZ2.open( params.read, false ) : FileCompression.NONE.open( params.read, false );
            map_stream  = params.bz2 ? FileCompression.BZ2.open( params.map, false ) : FileCompression.NONE.open( params.map, false );

            if( !te.read( read_stream, map_stream ) )
                throw new Exception( "fail" );
            
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
