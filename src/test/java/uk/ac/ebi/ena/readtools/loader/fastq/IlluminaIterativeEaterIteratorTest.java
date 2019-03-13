package uk.ac.ebi.ena.readtools.loader.fastq;

import java.io.File;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import uk.ac.ebi.ena.readtools.loader.common.QualityNormalizer;
import uk.ac.ebi.ena.readtools.loader.fastq.IlluminaIterativeEater.READ_TYPE;



public class 
IlluminaIterativeEaterIteratorTest
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
    
    //run ERR246213
    @Test
    public void 
    iteratorPairedTest1()
    {
        IlluminaIterativeEater wrapper = new IlluminaIterativeEater();
        
        wrapper.setFiles( new File[] { new File( "resources/T966_R1.fastq.gz" ),  
                                       new File( "resources/T966_R2.fastq.gz" ) } );
        
        wrapper.setNormalizers( new QualityNormalizer[] { QualityNormalizer.SANGER } );
        wrapper.setReadType( READ_TYPE.PAIRED );
       
        
        int spot_count = 0;
        int base_count = 0;
        
        Map<Integer, String> picked_reads = new HashMap<Integer, String>(); 
               
        //[era@ebi-025 ~]$ vdb-dump -R 1,759,2064,3829 -C READ /nfs/era-pub/vol1/err/ERR246/ERR246213
        picked_reads.put( 1,    "CATACACACTATAACAATAATGTCTATACTCACTAATTTTAGAATAAAACTTTAAACATTTATCACATACAGCATATGGATTCCCATCTCTATATACTATGCCAGAAAGTTACCACAGTTATGCACAGAGCTGCAAACAACTATACATGATATAATATTAGAATGTGTGTACTGCAAGC" );
        picked_reads.put( 759,  "CATACACACTATAACAATAATGTCTATACTCACTAATTTTAGAATAAAACTTTAAACATTTATCACATACAGCATATGGATTCCCATCTCTATATACTATGCCAGAAAGTTACCACAGTTATGCACAGAGCTGCAAACAACTATACATGATATAATATTAGAATGTGTGTACTGCAAGC" );
        picked_reads.put( 2064, "CATACACACTATAACAATAATGTCTATACTCACTAATTTTAGAATAAAACTTTAAACATTTATCACATACAGCATATGGATTCCCATCTCTATATACTATGCCAGAAAGTTACCACAGTTATGCACAGAGCTGCAAACAACTATACATGATATAATATTAGAATGTGTGTACTGCAAGC" );
        picked_reads.put( 3829, "CATACACACTATAACAATAATGTCTATACTCACTAATTTTAGAATAAAACTTTAAACATTTATCACATACAGCATATGGATTCCCATCTCTATATACTATGCCAGAAAGTTACCACAGTTATGCACAGAGCTGCAAACAACTATACATGATATAATATTAGAATGTGTGTACTGCAAGC" );

        
        for( Iterator<IlluminaSpot> i = wrapper.iterator(); i.hasNext(); )
        {
            IlluminaSpot is = i.next();
            base_count += is.bases.length();
            spot_count += 1;
            if( picked_reads.containsKey( spot_count ) )
                if( !picked_reads.get( spot_count ).equals( is.bases ) )
                    throw new RuntimeException( "Read " + spot_count + " failed" );
        }
        
        System.out.printf( "read_count: %d, base_count: %d\n", spot_count, base_count );
        if( 3829 != spot_count || base_count != 685391 )
            throw new RuntimeException( "" );
        System.out.println( "passed" );
    }


    @Test
    public void 
    iteratorPairedTest2()
    {
        IlluminaIterativeEater wrapper = new IlluminaIterativeEater();
        
        wrapper.setFiles( new File[] { new File( "bin/uk/ac/ebi/ena/frankenstein/loader/fastq/fastq_spots_correct2_1.txt" ),  
                                       new File( "bin/uk/ac/ebi/ena/frankenstein/loader/fastq/fastq_spots_correct2_2.txt" ) } );
        
        wrapper.setNormalizers( new QualityNormalizer[] { QualityNormalizer.SANGER } );
        wrapper.setReadType( READ_TYPE.PAIRED );
        
        int spot_count = 0;
        int base_count = 0;
        
        for( Iterator<IlluminaSpot> i = wrapper.iterator(); i.hasNext(); )
        {
            IlluminaSpot is = i.next();
            base_count += is.bases.length();
            spot_count += 1;
        }
        
        System.out.printf( "read_count: %d, base_count: %d\n", spot_count, base_count );
        if( 3 != spot_count || base_count != 606 )
            throw new RuntimeException( "" );
        System.out.printf( "passed\n" );
    }

}
