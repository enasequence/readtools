/*
* Copyright 2010-2021 EMBL - European Bioinformatics Institute
* Licensed under the Apache License, Version 2.0 (the "License"); you may not use this
* file except in compliance with the License. You may obtain a copy of the License at
* http://www.apache.org/licenses/LICENSE-2.0
* Unless required by applicable law or agreed to in writing, software distributed under the
* License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
* CONDITIONS OF ANY KIND, either express or implied. See the License for the
* specific language governing permissions and limitations under the License.
*/
package uk.ac.ebi.ena.readtools.loader.fastq;

import java.io.File;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

import org.junit.After;
import org.junit.Before;
import org.junit.Ignore;
import org.junit.Test;

import uk.ac.ebi.ena.readtools.common.reads.QualityNormalizer;
import uk.ac.ebi.ena.readtools.common.reads.normalizers.htsjdk.StandardQualityNormalizer;
import uk.ac.ebi.ena.readtools.loader.fastq.FastqIterativeWriter.READ_TYPE;



public class
MultiFastqConverterIteratorTest
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
        FastqIterativeWriter wrapper = new FastqIterativeWriter();
        
        wrapper.setFiles( new File[] { new File( "resources/T966_R1.fastq.gz" ),  
                                       new File( "resources/T966_R2.fastq.gz" ) } );
        
        wrapper.setNormalizers( new QualityNormalizer[] { new StandardQualityNormalizer() } );
        wrapper.setReadType( READ_TYPE.PAIRED );
       
        
        int spot_count = 0;
        int base_count = 0;
        
        Map<Integer, String> picked_reads = new HashMap<Integer, String>(); 
               
        //[era@ebi-025 ~]$ vdb-dump -R 1,759,2064,3829 -C READ /nfs/era-pub/vol1/err/ERR246/ERR246213
        picked_reads.put( 1,    "CATACACACTATAACAATAATGTCTATACTCACTAATTTTAGAATAAAACTTTAAACATTTATCACATACAGCATATGGATTCCCATCTCTATATACTATGCCAGAAAGTTACCACAGTTATGCACAGAGCTGCAAACAACTATACATGATATAATATTAGAATGTGTGTACTGCAAGC" );
        picked_reads.put( 759,  "CATACACACTATAACAATAATGTCTATACTCACTAATTTTAGAATAAAACTTTAAACATTTATCACATACAGCATATGGATTCCCATCTCTATATACTATGCCAGAAAGTTACCACAGTTATGCACAGAGCTGCAAACAACTATACATGATATAATATTAGAATGTGTGTACTGCAAGC" );
        picked_reads.put( 2064, "CATACACACTATAACAATAATGTCTATACTCACTAATTTTAGAATAAAACTTTAAACATTTATCACATACAGCATATGGATTCCCATCTCTATATACTATGCCAGAAAGTTACCACAGTTATGCACAGAGCTGCAAACAACTATACATGATATAATATTAGAATGTGTGTACTGCAAGC" );
        picked_reads.put( 3829, "CATACACACTATAACAATAATGTCTATACTCACTAATTTTAGAATAAAACTTTAAACATTTATCACATACAGCATATGGATTCCCATCTCTATATACTATGCCAGAAAGTTACCACAGTTATGCACAGAGCTGCAAACAACTATACATGATATAATATTAGAATGTGTGTACTGCAAGC" );

        
        for(Iterator<PairedRead> i = wrapper.iterator(); i.hasNext(); )
        {
            PairedRead is = i.next();
//            base_count += is.bases.length();
            base_count += (is.forward.bases.length() + is.reverse.bases.length());
            spot_count += 1;
            if( picked_reads.containsKey( spot_count ) )
                if( !picked_reads.get( spot_count ).equals( is.forward.bases + is.reverse.bases ) )
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
        FastqIterativeWriter wrapper = new FastqIterativeWriter();
        
        wrapper.setFiles( new File[] { new File( "resources/fastq_spots_correct2_1a.txt" ),
                                       new File( "resources/fastq_spots_correct2_2a.txt" ) } );
        
        wrapper.setNormalizers( new QualityNormalizer[] { new StandardQualityNormalizer() } );
        wrapper.setReadType( READ_TYPE.PAIRED );
        
        int spot_count = 0;
        int base_count = 0;
        
        for(Iterator<PairedRead> i = wrapper.iterator(); i.hasNext(); )
        {
            PairedRead is = i.next();
            base_count += (is.forward.bases.length() + is.reverse.bases.length());
            spot_count += 1;
        }
        
        System.out.printf( "read_count: %d, base_count: %d\n", spot_count, base_count );
        if( 3 != spot_count || base_count != 606 )
            throw new RuntimeException( "" );
        System.out.printf( "passed\n" );
    }

    @Test
    public void
    iteratorPairedTest3()
    {
        FastqIterativeWriter wrapper = new FastqIterativeWriter();

        wrapper.setFiles( new File[] { new File( "resources/fastq_spots_correct2_1.txt" ),
                                       new File( "resources/fastq_spots_correct2_2.txt" ) } );

        wrapper.setNormalizers( new QualityNormalizer[] { new StandardQualityNormalizer() } );
        wrapper.setReadType( READ_TYPE.PAIRED );

        int spot_count = 0;
        int base_count = 0;

        for(Iterator<PairedRead> i = wrapper.iterator(); i.hasNext(); )
        {
            PairedRead is = i.next();
            base_count += (is.forward.bases.length() + is.reverse.bases.length());
            spot_count += 1;
        }

        System.out.printf( "read_count: %d, base_count: %d\n", spot_count, base_count );
        if( 3 != spot_count || base_count != 606 )
            throw new RuntimeException( "" );
        System.out.printf( "passed\n" );
    }

}
