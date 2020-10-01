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
package uk.ac.ebi.ena.readtools.loader.bam;

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

import uk.ac.ebi.ena.readtools.cram.ref.ENAReferenceSource;
import uk.ac.ebi.ena.readtools.loader.common.eater.DataConsumer;
import uk.ac.ebi.ena.readtools.loader.common.eater.DataConsumerException;
import uk.ac.ebi.ena.readtools.loader.common.eater.PrintDataEater;
import uk.ac.ebi.ena.readtools.loader.common.feeder.DataProducer;
import uk.ac.ebi.ena.readtools.loader.common.feeder.DataProducerEOFException;
import uk.ac.ebi.ena.readtools.loader.common.feeder.DataProducerPanicException;

public class 
BamFeeder extends Thread implements DataProducer<BamSpot>
{
    SamReader           reader = null;
    Iterator<SAMRecord> it     = null;
    DataConsumer<BamSpot, ?> eater = new PrintDataEater<BamSpot, Void>();
    volatile boolean    is_ok = true;
    Throwable           stored_exception;
    
    
    public 
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
    produce() throws DataProducerPanicException, DataProducerEOFException
    {
        while( it.hasNext() ) 
            return new BamSpot( it.next() );
        try
        {
            reader.close();
        } catch( IOException ioe )
        {
            throw new DataProducerPanicException( ioe );
        }
        
        throw new DataProducerEOFException();
    }

    
    @Override
    public DataProducer<BamSpot>
    setConsumer(DataConsumer<BamSpot, ?> consumer)
    {
        this.eater = consumer;
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
                record = produce();
                ++ i;
                eater.consume( record );
            }while( true );
            
        } catch( DataProducerEOFException eof )
        {
            System.out.println( "EOF. Records: " + i );
            
        } catch( DataConsumerException | IllegalArgumentException e )
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
