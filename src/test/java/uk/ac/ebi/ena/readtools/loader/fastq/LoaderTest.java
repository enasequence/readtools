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
import java.io.InputStream;
import java.nio.file.Files;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import uk.ac.ebi.ena.readtools.loader.common.converter.AutoNormalizeQualityReadConverter;
import uk.ac.ebi.ena.readtools.loader.common.converter.ConverterException;

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
    read(File file, InputStream is) throws SecurityException, ConverterException {
        AutoNormalizeQualityReadConverter converter
                = new AutoNormalizeQualityReadConverter(is, new PrintReadWriter<>(), "", file.toString());
        try {
            converter.run();
            return true;
        } catch (Exception e) {
            e.printStackTrace();
            return false;
        }
    }
    
    
    boolean 
    read(String resource) throws Exception {
        File file = new File("src/test/resources/" + resource);
        InputStream is = getClass().getResourceAsStream( "src/test/resources/" + resource );
        is = is == null ? Files.newInputStream(file.toPath()) : is;
        try {
            return read(file, is);
        } finally {
            is.close();
        }
    }

    
    @Test public void
    testCorrect() throws Exception
    {
        if( !read("mp3_schw3.fq") )
            throw new Exception( "fail!" );

        if( !read("fastq_spots_correct.txt" ) )
            throw new Exception( "fail!" );

        if( !read("fastq_casava1_8_correct.txt" ) )
            throw new Exception( "fail!" );

    }
    
    
    @Test public void
    testFailed() throws Exception
    {
        if( read("fastq_spot_incorrect.txt") )
            throw new Exception( "fail!" );

        if( read("fastq_spot_incorrect2.txt") )
            throw new Exception( "fail!" );

        if( read("fastq_spot_incorrect3.txt") )
            throw new Exception( "fail!" );
        
        if( read("fastq_spot_incorrect4.txt") )
            throw new Exception( "fail!" );

        if( read("fastq_spot_incorrect5.txt") )
            throw new Exception( "fail!" );

        if( read("fastq_spot_incorrect6.txt") )
            throw new Exception( "fail!" );

        if( read("fastq_casava1_8_incorrect.txt") )
            throw new Exception( "fail!" );

    }
    
    
    
}
