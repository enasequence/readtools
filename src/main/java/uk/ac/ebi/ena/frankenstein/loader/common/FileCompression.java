package uk.ac.ebi.ena.frankenstein.loader.common;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.util.zip.GZIPInputStream;
import java.util.zip.ZipInputStream;

import org.apache.commons.compress.archivers.tar.TarArchiveInputStream;
import org.apache.commons.compress.compressors.bzip2.BZip2CompressorInputStream;


public enum 
FileCompression
{
    BZ2,
    GZIP,
    GZ,
    ZIP,
    BGZIP,
    BGZ,
    NONE;
    

    public InputStream
    open( String f_name, boolean tar ) throws FileNotFoundException, IOException
    {
        return open( new File( f_name ), tar );
    }
    

    public InputStream
    open( File file, boolean tar ) throws FileNotFoundException, IOException
    {
        InputStream is = new FileInputStream( file );
        try
        {
            switch( this )
            {
            case BZ2:
                is = new BZip2CompressorInputStream( new BufferedInputStream( is ), true );
                break;
            case GZIP:
            case GZ:
                is =  new GZIPInputStream( new BufferedInputStream( is ) );
                break;
            case ZIP:
                is = new ZipInputStream( new BufferedInputStream( is ) );
                break;
    
            case NONE:
                break;
            default:
                is.close();
                throw new RuntimeException( "Not implemented" );
            }
            
            if( tar ) 
            {
                @SuppressWarnings("resource")
                TarArchiveInputStream tais = new TarArchiveInputStream( is );
                tais.getNextTarEntry();
                is = tais;
            }
            return new BufferedInputStream( is );
        } catch( IOException e )
        {
            try
            {
                is.close();
            } catch( IOException e1 )
            {
                e1.printStackTrace();
            }
            
            throw e;
        }
    }

    
    public static FileCompression
    getCompressor( File file )
    {
        for( FileCompression c : FileCompression.values() )
            if( file.getName().toUpperCase().endsWith( "." + c.toString() ) )
                return c;
        return FileCompression.NONE;
    }
    
    
    public static InputStream
    open( File file ) throws FileNotFoundException, IOException
    {
        return getCompressor( file ).open( file.getPath(), getUseTar( file ) );
    }
    
    
    public static boolean
    getUseTar( File file )
    {
        return file.getName().toLowerCase().contains( ".tar" );
    }
}
