package uk.ac.ebi.ena.frankenstein.loader.bam;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;

import org.junit.After;
import org.junit.Before;

import uk.ac.ebi.ena.frankenstein.loader.common.eater.DataEaterException;
import uk.ac.ebi.ena.frankenstein.loader.common.eater.PrintDataEater;
import uk.ac.ebi.ena.frankenstein.loader.common.feeder.DataFeederException;
import uk.ac.ebi.ena.frankenstein.loader.fastq.IlluminaSpot;

public class 
LoaderTest
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
    read( InputStream is, String name ) throws SecurityException, DataFeederException, NoSuchMethodException, IOException, InterruptedException
    {
        BamEater eater = new BamEater( new File( "." ), 
                                       false, //read_type == IlluminaVDB2Eater.READ_TYPE.PAIRED, 
                                       true, 
                                       2 );
        eater.setVerbose( true );
        eater.setEater( new PrintDataEater<IlluminaSpot, Object>() );
        BamFeeder feeder = new BamFeeder( is ); 
        feeder.setName( name );
        feeder.setEater( eater );
        feeder.start();
        
        feeder.join();
        return feeder.isOk();
    }
    
    
    boolean 
    read( InputStream is1, 
          InputStream is2, 
          String name ) throws SecurityException, DataFeederException, NoSuchMethodException, IOException, InterruptedException, DataEaterException
    {
        BamEater eater = new BamEater( new File( "." ), 
                                       true, 
                                       true, 
                                       2 );
        //eater.setVerbose( true );
        eater.setEater( new PrintDataEater<IlluminaSpot, Object>() );
        
        BamFeeder feeder1 = new BamFeeder( is1 ); 
        feeder1.setName( name + ".1" );
        feeder1.setEater( eater );
        feeder1.start();
        /*
        BamFeeder feeder2 = new BamFeeder( is2 ); 
        feeder2.setName( name + ".2" );
        feeder2.setEater( eater );
        feeder2.start();
         */
        feeder1.join();
        //feeder2.join();
        
        eater.cascadeErrors();
        return eater.isOk() && feeder1.isOk();// && feeder2.isOk();
    }

    
    boolean 
    read( File resource1, 
          File resource2,
          String name ) throws Exception
    {
        InputStream is1 = new BufferedInputStream( new FileInputStream( resource1 ) );
        InputStream is2 = new BufferedInputStream( new FileInputStream( resource2 ) );
        try
        {
            return read( is1, is2, name );
        } finally
        {
            is1.close();
            is2.close();
        }
    }
    
    
    boolean 
    read( String resource ) throws Exception
    {
        InputStream is = getClass().getResourceAsStream( resource );
        try
        {
            return read( is, resource );
        } finally
        {
            is.close();
        }
    }
    
    
    boolean 
    read( File file ) throws Exception
    {
        InputStream is = new BufferedInputStream( new FileInputStream( file ) );
        try
        {
            return read( is, file.getPath() );
        } finally
        {
            is.close();
        }
    }

    
    
    @org.junit.Test
    public void
    testCorrect() throws Exception
    {
        if( !read( new File( "src/test/java/uk/ac/ebi/ena/frankenstein/loader/bam/bnlx1mp_srt.header.bam" ),
                   new File( "src/test/java/uk/ac/ebi/ena/frankenstein/loader/bam/bnlx1mp_srt.header.bam" ), 
                   "bnlx1mp_srt.header.bam" ) )
            throw new Exception( "fail!" );
       
    }
    
    
    @org.junit.Test
    public void
    testFailed() throws Exception
    {
//        if( read( "fastq_casava1_8_incorrect.txt", QualityNormalizer.SANGER ) )
//            throw new Exception( "fail!" );

    }
    
    
    
}
