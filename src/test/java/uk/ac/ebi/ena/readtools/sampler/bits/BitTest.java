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
package uk.ac.ebi.ena.readtools.sampler.bits;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.EOFException;
import java.io.IOException;

import org.junit.Assert;
import org.junit.Test;



public class BitTest
{
    //@Test
    public void
    BitStreamTest() throws IOException
    {
        ByteArrayOutputStream bytes = new ByteArrayOutputStream();
        BitOutputStream bstream = new BitOutputStream( bytes );
        bstream.write( 0, 3 );
        bstream.write( 1, 3 );
        bstream.write( 2, 3 );
        bstream.write( 3, 3 );
        bstream.write( 4, 3 );
        bstream.write( 5, 3 );
        bstream.write( 6, 3 );
        bstream.write( 7, 3 );
        bstream.write( 265535, 21 );
        bstream.write( 265536, 23 );
        bstream.write(  1073741836, 31 );
        bstream.close();
        
        BitInputStream bis = new BitInputStream( new ByteArrayInputStream( bytes.toByteArray() ) );
        System.out.println( "Stream length: " + bytes.toByteArray().length );
        System.out.printf( "%d, %d, %d, %d, %d, %d, %d, %d\n", 
                           bis.read( 3 ), bis.read( 3 ), bis.read( 3 ), bis.read( 3 ),
                           bis.read( 3 ), bis.read( 3 ), bis.read( 3 ), bis.read( 3 ) );
        System.out.println( bis.read( 21 ) );
        
        System.out.println( bis.read( 23 ) );
        System.out.println( bis.read( 31 ) );
        
        bis.close();
        
    }


    @Test
    public void
    DictStreamTest() throws IOException
    {
        ByteArrayOutputStream bytes = new ByteArrayOutputStream();
        DictionaryOutputStream bstream = new DictionaryOutputStream( new BitOutputStream( bytes ), new BaseDict() );
        String seq = "ACGTNTGCAN";
        bstream.write( seq.getBytes() );
        
        bstream.close();
        
        DictionaryInputStream bis = new DictionaryInputStream( new BitInputStream( new ByteArrayInputStream( bytes.toByteArray() ) ), new BaseDict() ); 
        System.out.println( "Stream length: " + bytes.toByteArray().length );

        StringBuilder sb = new StringBuilder();
        int ch = 0;
        
        try
        {
        	while( -1 != ( ch = bis.read() ) )
        		sb.append( (char)ch );

        	Assert.fail();
        } catch( EOFException eof )
        {
        	Assert.assertEquals( seq, sb.toString() );
        }
    }
}
