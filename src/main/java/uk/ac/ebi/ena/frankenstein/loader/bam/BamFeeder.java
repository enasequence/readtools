package uk.ac.ebi.ena.frankenstein.loader.bam;

import java.io.IOException;
import java.io.InputStream;
import java.util.Iterator;

import htsjdk.samtools.DefaultSAMRecordFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.Log.LogLevel;
import net.sf.cram.ref.ENAReferenceSource;
import uk.ac.ebi.ena.frankenstein.loader.common.eater.DataEater;
import uk.ac.ebi.ena.frankenstein.loader.common.eater.DataEaterException;
import uk.ac.ebi.ena.frankenstein.loader.common.eater.PrintDataEater;
import uk.ac.ebi.ena.frankenstein.loader.common.feeder.DataFeeder;
import uk.ac.ebi.ena.frankenstein.loader.common.feeder.DataFeederEOFException;
import uk.ac.ebi.ena.frankenstein.loader.common.feeder.DataFeederPanicException;

public class 
BamFeeder extends Thread implements DataFeeder<BamSpot>
{
    SamReader           reader = null;
    Iterator<SAMRecord> it     = null;
    DataEater<BamSpot, ?> eater = new PrintDataEater<BamSpot, Void>();
    volatile boolean    is_ok = true;
    Throwable           stored_exception;
    
    
    protected 
    BamFeeder( InputStream istream ) 
    {
        Log.setGlobalLogLevel( LogLevel.ERROR );
        SamReaderFactory.setDefaultValidationStringency( ValidationStringency.SILENT );
        SamReaderFactory factory = SamReaderFactory.make();
        factory.enable( SamReaderFactory.Option.DONT_MEMORY_MAP_INDEX );
        factory.validationStringency( ValidationStringency.SILENT );
        factory.referenceSource( new ENAReferenceSource() );
        factory.samRecordFactory( DefaultSAMRecordFactory.getInstance() );
        SamInputResource ir = SamInputResource.of( istream );

        //File indexMaybe = SamFiles.findIndex(file);
        //reporter.write(Severity.INFO, "proposed index: " + indexMaybe);

//        if (null != indexMaybe)
//            ir.index(indexMaybe);

        reader = factory.open( ir );

        System.out.printf( "BAM header version %s\n", reader.getFileHeader().getAttributes() );
        it = reader.iterator();
    }


    @Override
    public BamSpot 
    feed() throws DataFeederPanicException, DataFeederEOFException 
    {
        while( it.hasNext() ) 
            return new BamSpot( it.next() );
        try
        {
            reader.close();
        } catch( IOException ioe )
        {
            throw new DataFeederPanicException( ioe );
        }
        
        throw new DataFeederEOFException();
    }

    
    @Override
    public DataFeeder<BamSpot> 
    setEater( DataEater<BamSpot, ?> eater )
    {
        this.eater = eater;
        return this;
    }

    
    public void
    run()
    {
        int i = 0;
        BamSpot record;
        try
        {
            do
            {
                record = feed();
                ++ i;
                eater.eat( record );
            }while( true );
            
        } catch( DataFeederEOFException eof )
        {
            System.out.println( "EOF. Records: " + i );
            
        } catch( DataEaterException | IllegalArgumentException e )
        {
            
            e.printStackTrace();
            is_ok = false;
            stored_exception = e;
            
        } catch( Throwable t )
        {
            System.out.println( "Failed on record: " + i );
            t.printStackTrace();
            is_ok = false;
            stored_exception = t;
        }
    }
    
    
    public boolean
    isOk()
    {
        return is_ok;
    }
    
    
    public Throwable
    getStoredException()
    {
        return stored_exception;
    }
    
}
