package uk.ac.ebi.ena.frankenstein.loader.fastq;

import java.io.File;
import java.util.Iterator;

import uk.ac.ebi.ena.frankenstein.loader.common.QualityNormalizer;

public class 
IlluminaIterativeEater implements Iterable<IlluminaSpot>
{
    public enum
    READ_TYPE
    {
      SINGLE,
      PAIRED
    };
    
    
    private File[] files;
    private int spill_page_size = 5000000;
    private File tmp_folder = new File( "." );
    private READ_TYPE read_type;
    private QualityNormalizer[] normalizers;
    
    
    @Override
    public Iterator<IlluminaSpot> 
    iterator()
    {
        try
        {
            return (Iterator<IlluminaSpot>) new IlluminaIterativeEaterIterator( tmp_folder, 
                                                                        spill_page_size, 
                                                                        read_type, 
                                                                        files,
                                                                        normalizers );
        } catch( Throwable t )
        {
            t.printStackTrace();
        }
        
        return null;
    }


    public void
    setFiles( File files[] )
    {
        this.files = files;
    }

    
    public File[]
    getFiles()
    {
        return this.files;
    }
    
    
    public int 
    getSpillPageSize()
    {
        return spill_page_size;
    }


    public void 
    setSpillPageSize( int spill_page_size )
    {
        this.spill_page_size = spill_page_size;
    }


    public File 
    getTmpFolder()
    {
        return tmp_folder;
    }


    public void 
    setTmpFolder( File tmp_folder )
    {
        this.tmp_folder = tmp_folder;
    }


    public READ_TYPE 
    getReadType()
    {
        return read_type;
    }


    public void 
    setReadType( READ_TYPE read_type )
    {
        this.read_type = read_type;
    }


    public QualityNormalizer[] 
    getNormalizers()
    {
        return normalizers;
    }


    public void 
    setNormalizers( QualityNormalizer[] normalizers )
    {
        this.normalizers = normalizers;
    }
}
