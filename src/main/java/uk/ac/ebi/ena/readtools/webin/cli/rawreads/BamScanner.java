package uk.ac.ebi.ena.readtools.webin.cli.rawreads;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.stream.Collectors;

import org.apache.log4j.Logger;

import htsjdk.samtools.DefaultSAMRecordFactory;
import htsjdk.samtools.SAMFormatException;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamFiles;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.cram.CRAMException;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.Log.LogLevel;
import uk.ac.ebi.ena.readtools.cram.ref.ENAReferenceSource;
import uk.ac.ebi.ena.readtools.webin.cli.rawreads.ScannerMessage.ScannerErrorMessage;
import uk.ac.ebi.ena.readtools.webin.cli.rawreads.ScannerMessage.ScannerInfoMessage;
import uk.ac.ebi.ena.readtools.webin.cli.rawreads.refs.CramReferenceInfo;


public abstract class
BamScanner 
{
    static 
    {
        System.setProperty( "samjdk.use_cram_ref_download", Boolean.TRUE.toString() );
    }


    private static final Logger log = Logger.getLogger( BamScanner.class );
    private static final String BAM_STAR = "*";    
    abstract protected void logProcessedReadNumber( long cnt );
    
    
    public List<ScannerMessage>
    readCramFile( RawReadsFile rf, AtomicBoolean paired ) throws IOException
    {
        List<ScannerMessage> result = new ArrayList<>();
        CramReferenceInfo cri = new CramReferenceInfo();
        {
            Map<String, Boolean> ref_set;
            try
            {
                ref_set = cri.confirmFileReferences( new File( rf.getFilename() ) );
                if( !ref_set.isEmpty() && ref_set.containsValue( Boolean.FALSE ) )
                {
                    result.add( new ScannerErrorMessage( "Unable to find reference sequence(s) from the CRAM reference registry: " 
                                                       + ref_set.entrySet()
                                                                .stream()
                                                                .filter( e -> !e.getValue() )
                                                                .map( e -> e.getKey() )
                                                                .collect(Collectors.toList() ) ) );
                }

            } catch( IOException e1 )
            {
                result.add( new ScannerErrorMessage( e1, e1.getMessage(), null ) );
            }
        }
        
        result.addAll( readBamFile( rf, paired ) );
        return result;
    }
    
    
    public List<ScannerMessage>
    readBamFile( RawReadsFile rf, AtomicBoolean paired ) throws IOException
    {
        List<ScannerMessage> reporter = new ArrayList<>();
        {
            long read_no = 0;
            long reads_cnt = 0;

            try 
            {
                String msg = String.format( "Processing file %s\n", rf.getFilename() );
                log.info( msg );
                
                ENAReferenceSource reference_source = new ENAReferenceSource();
                reference_source.setLoggerWrapper( new ENAReferenceSource.LoggerWrapper() 
                {
                    @Override public void
                    error( Object... messageParts )
                    {
                        reporter.add( new ScannerErrorMessage( null == messageParts ? "null" : String.valueOf( Arrays.asList( messageParts ) ) ) );
                    }

                    @Override public void
                    warn( Object... messageParts )
                    {
                        reporter.add( new ScannerInfoMessage( /* Severity.WARNING, */ null == messageParts ? "null" : String.valueOf( Arrays.asList( messageParts ) ) ) );
                    }

                    @Override public void
                    info( Object... messageParts )
                    {
                        reporter.add( new ScannerInfoMessage( null == messageParts ? "null" : String.valueOf( Arrays.asList( messageParts ) ) ) );
                    }
                } );

                reporter.add( new ScannerInfoMessage( "REF_PATH  " + reference_source.getRefPathList() ) );
                reporter.add( new ScannerInfoMessage( "REF_CACHE " + reference_source.getRefCacheList() ) );


                File file = new File( rf.getFilename() );
                Log.setGlobalLogLevel( LogLevel.ERROR );
                SamReaderFactory.setDefaultValidationStringency( ValidationStringency.SILENT );
                SamReaderFactory factory = SamReaderFactory.make();
                factory.enable( SamReaderFactory.Option.DONT_MEMORY_MAP_INDEX );
                factory.validationStringency( ValidationStringency.SILENT );
                factory.referenceSource( reference_source );
                factory.samRecordFactory( DefaultSAMRecordFactory.getInstance() );
                SamInputResource ir = SamInputResource.of( file );
                File indexMaybe = SamFiles.findIndex( file );
                reporter.add( new ScannerInfoMessage( "proposed index: " + indexMaybe ) );

                if (null != indexMaybe)
                    ir.index(indexMaybe);

                SamReader reader = factory.open(ir);

                for( SAMRecord record : reader )
                {
                    read_no++;
                    //do not load supplementary reads
                    if( record.isSecondaryOrSupplementary() )
                        continue;

                    if( record.getDuplicateReadFlag() )
                        continue;

                    if( record.getReadString().equals( BAM_STAR ) && record.getBaseQualityString().equals( BAM_STAR ) )
                        continue;

                    if( record.getReadBases().length != record.getBaseQualities().length )
                    {
                        ScannerMessage validationMessage = new ScannerErrorMessage( null, 
                                                                                    "Mismatch between length of read bases and qualities",
                                                                                    String.format( "%s:%d", rf.getFilename(), read_no ) );
                        
                        reporter.add( validationMessage );
                    }

                    paired.compareAndSet( false, record.getReadPairedFlag() );
                    reads_cnt++;
                    if( 0 == reads_cnt % 1000 )
                        logProcessedReadNumber( reads_cnt );
                }

                logProcessedReadNumber( reads_cnt );

                reader.close();

                reporter.add( new ScannerInfoMessage( "Valid reads count: " + reads_cnt ) );
                reporter.add( new ScannerInfoMessage( "LibraryLayout: " + ( paired.get() ? "PAIRED" : "SINGLE" ) ) );

                if( 0 == reads_cnt )
                    reporter.add( new ScannerErrorMessage( "File contains no valid reads" ) );

            } catch( SAMFormatException | CRAMException e )
            {
                reporter.add( new ScannerErrorMessage( e.getMessage() ) );

            }
        }
        return reporter;
    }
}
