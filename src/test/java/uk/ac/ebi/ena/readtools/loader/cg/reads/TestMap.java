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
        InputStream is = TestMap.class.getClassLoader().getResourceAsStream( resource );
        
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
