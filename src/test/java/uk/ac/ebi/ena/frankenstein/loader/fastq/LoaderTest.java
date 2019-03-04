package uk.ac.ebi.ena.frankenstein.loader.fastq;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.concurrent.atomic.AtomicLong;
import java.util.concurrent.atomic.AtomicReference;

import org.junit.After;
import org.junit.Before;

import uk.ac.ebi.ena.frankenstein.loader.common.QualityNormalizer;
import uk.ac.ebi.ena.frankenstein.loader.common.eater.PrintDataEater;
import uk.ac.ebi.ena.frankenstein.loader.common.feeder.AbstractDataFeeder;
import uk.ac.ebi.ena.frankenstein.loader.common.feeder.DataFeederException;
import uk.ac.ebi.ena.frankenstein.loader.fastq.DataSpot.DataSpotParams;

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
    read( InputStream is, String name, final QualityNormalizer normalizer ) throws SecurityException, DataFeederException, NoSuchMethodException, IOException, InterruptedException
    {
        AbstractDataFeeder<DataSpot> df = new AbstractDataFeeder<DataSpot>( is, DataSpot.class ) 
        {
            final AtomicLong line_no = new AtomicLong( 1 );
            final AtomicReference<DataSpot.ReadStyle> read_style = new AtomicReference<DataSpot.ReadStyle>();
            DataSpotParams params = DataSpot.defaultParams();
            @Override
            protected DataSpot newFeedable()
            {
                return new DataSpot( normalizer, "", params );
            }
        };
        df.setName( name );
        df.setEater( new PrintDataEater<DataSpot, Object>() );
        df.start();
        df.join();
        return df.isOk();
    }
    
    
    boolean 
    read( String resource, final QualityNormalizer normalizer ) throws Exception
    {
        InputStream is = getClass().getResourceAsStream( resource );
        try
        {
            return read( is, resource, normalizer );
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
            return read( is, file.getPath(), normalizer );
        } finally
        {
            is.close();
        }
    }

    
    
    @org.junit.Test
    public void
    testCorrect() throws Exception
    {
        
        if( !read( new File( "src/test/java/uk/ac/ebi/ena/frankenstein/loader/fastq/mp3_schw3.fq" ), QualityNormalizer.ILLUMINA_1_3 ) )
            throw new Exception( "fail!" );
        
        if( !read( new File( "src/test/java/uk/ac/ebi/ena/frankenstein/loader/fastq/fastq_spots_correct.txt" ), QualityNormalizer.ILLUMINA_1_3 ) )
            throw new Exception( "fail!" );
        
        if( !read( new File( "src/test/java/uk/ac/ebi/ena/frankenstein/loader/fastq/fastq_casava1_8_correct.txt" ), QualityNormalizer.SANGER ) )
            throw new Exception( "fail!" );

    }
    
    
    @org.junit.Test
    public void
    testFailed() throws Exception
    {

        if( read( new File( "src/test/java/uk/ac/ebi/ena/frankenstein/loader/fastq/fastq_spot_incorrect.txt" ), QualityNormalizer.ILLUMINA_1_3 ) )
            throw new Exception( "fail!" );

        if( read( new File( "src/test/java/uk/ac/ebi/ena/frankenstein/loader/fastq/fastq_spot_incorrect2.txt" ), QualityNormalizer.ILLUMINA_1_3 ) )
            throw new Exception( "fail!" );

        if( read( new File( "src/test/java/uk/ac/ebi/ena/frankenstein/loader/fastq/fastq_spot_incorrect3.txt" ), QualityNormalizer.ILLUMINA_1_3 ) )
            throw new Exception( "fail!" );
        
        if( read( new File( "src/test/java/uk/ac/ebi/ena/frankenstein/loader/fastq/fastq_spot_incorrect4.txt" ), QualityNormalizer.ILLUMINA_1_3 ) )
            throw new Exception( "fail!" );

        if( read( new File( "src/test/java/uk/ac/ebi/ena/frankenstein/loader/fastq/fastq_spot_incorrect5.txt" ), QualityNormalizer.ILLUMINA_1_3 ) )
            throw new Exception( "fail!" );
        
        if( read( new File( "src/test/java/uk/ac/ebi/ena/frankenstein/loader/fastq/fastq_spot_incorrect6.txt" ), QualityNormalizer.SANGER ) )
            throw new Exception( "fail!" );

        if( read( new File( "src/test/java/uk/ac/ebi/ena/frankenstein/loader/fastq/fastq_casava1_8_incorrect.txt" ), QualityNormalizer.SANGER ) )
            throw new Exception( "fail!" );

    }
    
    
    
}
