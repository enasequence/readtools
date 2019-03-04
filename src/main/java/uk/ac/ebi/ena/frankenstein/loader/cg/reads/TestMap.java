package uk.ac.ebi.ena.frankenstein.loader.cg.reads;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;

import org.junit.After;
import org.junit.Before;

import uk.ac.ebi.ena.frankenstein.loader.common.QualityNormalizer;
import uk.ac.ebi.ena.frankenstein.loader.common.eater.PrintDataEater;
import uk.ac.ebi.ena.frankenstein.loader.common.feeder.AbstractDataFeeder;
import uk.ac.ebi.ena.frankenstein.loader.common.feeder.DataFeederException;

public class 
TestMap
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
        AbstractDataFeeder<CompleteGenomicsMap3> df = new AbstractDataFeeder<CompleteGenomicsMap3>( is, CompleteGenomicsMap3.class ) 
        {

            @Override
            protected CompleteGenomicsMap3 newFeedable()
            {
                return new CompleteGenomicsMap3( h );
            }
        };
        
        df.setEater( new PrintDataEater<CompleteGenomicsMap3, Object>() );
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
        if( !read( "cg_map.txt" ) )
            throw new Exception( "fail!" );
    }
    
    
}
