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

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;

import uk.ac.ebi.ena.readtools.loader.common.FileCompression;
import uk.ac.ebi.ena.readtools.loader.common.eater.DataConsumer;
import uk.ac.ebi.ena.readtools.loader.common.eater.DataConsumerException;
import uk.ac.ebi.ena.readtools.loader.common.eater.PrintDataEater;
import uk.ac.ebi.ena.readtools.loader.common.feeder.AbstractDataProducer;
import uk.ac.ebi.ena.readtools.loader.common.feeder.DataProducerException;

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
          InputStream map_stream ) throws SecurityException, DataProducerException, NoSuchMethodException, IOException, DataConsumerException
    {
        final ReadHeader read_header = ReadHeader( read_stream );
        final ReadHeader map_header  = ReadHeader( map_stream );
        if( !read_header.getFileProrerty( ReadHeader.PropertyName.BATCH_FILE_NUMBER ).equals( map_header.getFileProrerty( ReadHeader.PropertyName.BATCH_FILE_NUMBER ) ) )
                throw new IOException( "Batches did not match!" );
        
        AbstractDataProducer<CompleteGenomicsMap3> mf = new AbstractDataProducer<CompleteGenomicsMap3>( map_stream, CompleteGenomicsMap3.class )
        {

            @Override
            protected CompleteGenomicsMap3
            newProducible()
            {
                return new CompleteGenomicsMap3( map_header );
            }
        };

        AbstractDataProducer<CompleteGenomicsRead> rf = new AbstractDataProducer<CompleteGenomicsRead>( read_stream, CompleteGenomicsRead.class )
        {

            @Override
            protected CompleteGenomicsRead
            newProducible()
            {
                return new CompleteGenomicsRead( read_header );
            }
        };
        

        
        DataConsumer<? extends CompleteGenomicsBase, ?> eater = new CompleteGenomicsPairedEater();
        eater.setConsumer( new PrintDataEater() );
        mf.setConsumer( (DataConsumer<CompleteGenomicsMap3, ?>) eater );
        rf.setConsumer( (DataConsumer<CompleteGenomicsRead, ?>) eater );
        
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
