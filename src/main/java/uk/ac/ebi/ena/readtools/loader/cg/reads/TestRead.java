package uk.ac.ebi.ena.readtools.loader.cg.reads;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;

import org.junit.After;
import org.junit.Before;

import uk.ac.ebi.ena.readtools.loader.common.QualityNormalizer;
import uk.ac.ebi.ena.readtools.loader.common.eater.PrintDataEater;
import uk.ac.ebi.ena.readtools.loader.common.feeder.AbstractDataFeeder;
import uk.ac.ebi.ena.readtools.loader.common.feeder.DataFeederException;

public class 
TestRead
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
    read( InputStream is ) throws SecurityException, DataFeederException, NoSuchMethodException, IOException
    {
        final ReadHeader h = ReadHeader( is );
        
        AbstractDataFeeder<CompleteGenomicsRead> df = new AbstractDataFeeder<CompleteGenomicsRead>( is, CompleteGenomicsRead.class ) 
        {

            @Override
            protected CompleteGenomicsRead newFeedable()
            {
                return new CompleteGenomicsRead( h );
            }
        };
        
        df.setEater( new PrintDataEater<CompleteGenomicsRead, Object>() );
        df.run();
        return df.isOk();
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
    read( String resource ) throws Exception
    {
        InputStream is = getClass().getResourceAsStream( resource );
        
        try
        {
            return read( is );
        } finally
        {
            is.close();
        }
    }
    
    
    boolean 
    read( File file, final QualityNormalizer normalizer ) throws Exception
    {
        InputStream is = new BufferedInputStream( new FileInputStream( file ) );
        try
        {
            ReadHeader( is );
            return read( is );
        } finally
        {
            is.close();
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
        if( !read( "cg_read.txt" ) )
            throw new Exception( "fail!" );
    }
    
    
}
