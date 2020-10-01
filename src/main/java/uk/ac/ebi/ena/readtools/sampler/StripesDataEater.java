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
package uk.ac.ebi.ena.readtools.sampler;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.nio.ByteBuffer;

import org.apache.commons.compress.compressors.bzip2.BZip2CompressorOutputStream;

import uk.ac.ebi.ena.readtools.loader.common.eater.DataConsumer;
import uk.ac.ebi.ena.readtools.loader.common.eater.DataConsumerException;
import uk.ac.ebi.ena.readtools.loader.fastq.DataSpot;
import uk.ac.ebi.ena.readtools.sampler.bits.BaseDict;
import uk.ac.ebi.ena.readtools.sampler.bits.BitOutputStream;
import uk.ac.ebi.ena.readtools.sampler.bits.DictionaryOutputStream;
import uk.ac.ebi.ena.readtools.sampler.bits.ScoreDict;

public class 
StripesDataEater<T1 extends DataSpot, T2> implements DataConsumer<T1, T2>
{
    final static ByteBuffer reads_buf = ByteBuffer.allocateDirect( 4096 );
    final static ByteBuffer scores_buf = ByteBuffer.allocateDirect( 4096 );
    
    final protected int lbound = 33;  //!
    final protected int hbound = 126; //~;
    
    //final protected long row[] = new long[ hbound - lbound ] ;
    long  read_count;
    boolean lock;
    final File    folder;
    final File reads;
    final File scores;
    final File names;
    final File offsets;

    OutputStream r_stream;
    OutputStream s_stream[];
    OutputStream n_stream;
    OutputStream o_stream;
    BaseSqueezer bsq = new BaseSqueezer();
    
    public
    StripesDataEater( File folder ) throws IOException
    {
        this.folder = folder;
        
        if( !folder.exists() )
            folder.mkdirs();
        
        this.reads    = new File( folder, "reads" );
        if( !reads.exists() )
            reads.createNewFile();

        this.scores   = new File( folder, "scores" );
        if( !scores.exists() )
            scores.mkdir();
        
        this.names    = new File( folder, "names" );
        if( !names.exists() )
            names.createNewFile();
        
        this.offsets  = new File( folder, "offsets" );
        if( !offsets.exists() )
            offsets.createNewFile();
        
    }
    

    @Override
    public void 
    cascadeErrors() throws DataConsumerException
    {
        try
        {        
            if( null != r_stream )
                r_stream.close(); 
            
            if( null != s_stream )
                for( int i = 0; i < s_stream.length; ++i )
                    s_stream[ i ].close();

            if( null != n_stream )
                n_stream.close();

            if( null != o_stream )
                o_stream.close();

        } catch( IOException e )
        {
            throw new DataConsumerException( e );
        }
    }

   
    public boolean
    setLock()
    {
        if( lock )
            return false;
        lock = true;
        return lock;
    }
    

    public void
    resetLock()
    {
        lock = false;
    }


    @Override
    public void
    consume(T1 object ) throws DataConsumerException
    {
        try
        {        
            if( null == r_stream )
                r_stream = new DictionaryOutputStream( new BitOutputStream( new BZip2CompressorOutputStream( new BufferedOutputStream( new FileOutputStream( reads ) ) ) ), new BaseDict() );
            
            if( null == s_stream )
            {
                s_stream = new OutputStream[ object.quals.length() ];
                for( int i = 0; i < s_stream.length; ++i )
                    s_stream[ i ] = new DictionaryOutputStream( new BitOutputStream( new BZip2CompressorOutputStream( new BufferedOutputStream( new FileOutputStream( new File( scores, "_" + i ) ) ) ) ), new ScoreDict() );
            }
            if( null == n_stream )
                n_stream = new BZip2CompressorOutputStream( new BufferedOutputStream( new FileOutputStream( names ) ) );

            if( null == o_stream )
                o_stream = new BZip2CompressorOutputStream( new BufferedOutputStream( new FileOutputStream( offsets ) ) );

            r_stream.write( ( object.bases ).getBytes() );
            for( int i = 0; i < s_stream.length; ++i )
                s_stream[ i ].write( object.quals.charAt( i ) );
            
            n_stream.write( object.bname.getBytes() );
            
            o_stream.write( object.bases.length() );
            o_stream.write( ( object.bases.length() >> 8 ) );
            
            o_stream.write( object.quals.length() );
            o_stream.write( ( object.quals.length() >> 8 ) );
            
            o_stream.write( object.bname.length() );
            o_stream.write( ( object.bname.length() >> 8 ) );

            
        } catch( IOException e )
        {
            throw new DataConsumerException( e );
        }
        
        ++read_count;
    }


    @Override
    public void
    setConsumer(DataConsumer<T2, ?> dataConsumer)
    {
        throw new UnsupportedOperationException();
    }
    
    
    public long 
    getFeedCount()
    {
        return read_count;
    }
    
    
    public byte
    getLBoundary()
    {
        return lbound;
    }
    
    
    public byte
    getHBoundary()
    {
        return hbound;
    }


    @Override
    public boolean 
    isOk()
    {
        return true;
    }

}
