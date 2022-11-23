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

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStream;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import uk.ac.ebi.ena.readtools.common.reads.QualityNormalizer;
import uk.ac.ebi.ena.readtools.common.reads.normalizers.htsjdk.IlluminaQualityNormalizer;
import uk.ac.ebi.ena.readtools.common.reads.normalizers.htsjdk.StandardQualityNormalizer;
import uk.ac.ebi.ena.readtools.loader.common.converter.ConverterException;
import uk.ac.ebi.ena.readtools.loader.common.converter.ReadConverter;
import uk.ac.ebi.ena.readtools.loader.common.writer.PrintReadWriter;

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
    read( InputStream is, String name, final QualityNormalizer normalizer ) throws SecurityException, ConverterException, InterruptedException
    {
        ReadConverter df = new ReadConverter( is, normalizer, "" );
        df.setWriter( new PrintReadWriter<>() );
        df.run();
        return df.isOk();
    }
    
    
    boolean 
    read( String resource, final QualityNormalizer normalizer ) throws Exception
    {
        InputStream is = getClass().getResourceAsStream( "/resources/" + resource );
        is = is == null ? new FileInputStream( new File( "resources/" + resource ) ) : is;
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

    
    
    @Test public void
    testCorrect() throws Exception
    {
        IlluminaQualityNormalizer normalizer = new IlluminaQualityNormalizer();
        
        if( !read( "mp3_schw3.fq", normalizer) )
            throw new Exception( "fail!" );
        
        if( !read( "fastq_spots_correct.txt", normalizer ) )
            throw new Exception( "fail!" );
        
        if( !read( "fastq_casava1_8_correct.txt", normalizer ) )
            throw new Exception( "fail!" );

    }
    
    
    @Test public void
    testFailed() throws Exception
    {
        StandardQualityNormalizer stdNormalizer = new StandardQualityNormalizer();
        IlluminaQualityNormalizer illNormalizer = new IlluminaQualityNormalizer();

        if( read( "fastq_spot_incorrect.txt", illNormalizer ) )
            throw new Exception( "fail!" );

        if( read( "fastq_spot_incorrect2.txt", illNormalizer ) )
            throw new Exception( "fail!" );

        if( read( "fastq_spot_incorrect3.txt", illNormalizer ) )
            throw new Exception( "fail!" );
        
        if( read( "fastq_spot_incorrect4.txt", illNormalizer ) )
            throw new Exception( "fail!" );

        if( read( "fastq_spot_incorrect5.txt", illNormalizer ) )
            throw new Exception( "fail!" );

        if( read( "fastq_spot_incorrect6.txt", stdNormalizer ) )
            throw new Exception( "fail!" );

        if( read( "fastq_casava1_8_incorrect.txt", stdNormalizer ) )
            throw new Exception( "fail!" );

    }
    
    
    
}
