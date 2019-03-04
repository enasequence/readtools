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
