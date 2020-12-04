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
package uk.ac.ebi.ena.readtools.webin.cli.rawreads;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Map;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.stream.Collectors;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

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
import uk.ac.ebi.ena.readtools.webin.cli.rawreads.refs.CramReferenceInfo;
import uk.ac.ebi.ena.webin.cli.validator.message.ValidationMessage;
import uk.ac.ebi.ena.webin.cli.validator.message.ValidationOrigin;
import uk.ac.ebi.ena.webin.cli.validator.message.ValidationResult;


public abstract class
BamScanner 
{
    static 
    {
        System.setProperty( "samjdk.use_cram_ref_download", Boolean.TRUE.toString() );
    }

    //TODO remove duplication
    private static final int    DEFAULT_PRINT_FREQ = 1_000;
    private static final String PRINT_FREQ_PROPERTY_NAME = "webincli.scanner.print.freq";
    private static final int    print_freq = Integer.valueOf( System.getProperty( PRINT_FREQ_PROPERTY_NAME, String.valueOf( DEFAULT_PRINT_FREQ ) ) );    

    private static final Logger log = LoggerFactory.getLogger( BamScanner.class );
    private static final String BAM_STAR = "*";    
    abstract protected void logProcessedReadNumber( long cnt );


    /**
     *
     * @param fileResult
     * @param rf
     * @param paired
     * @param confirmCramReferences - Confirm CRAM file reference if set to true.
     * @throws IOException
     */
    public void
    readCramFile( ValidationResult fileResult, RawReadsFile rf, AtomicBoolean paired, boolean confirmCramReferences ) throws IOException
    {
        if (confirmCramReferences) {
            CramReferenceInfo cri = new CramReferenceInfo();
            {
                Map<String, Boolean> ref_set;
                try
                {
                    ref_set = cri.confirmFileReferences( new File( rf.getFilename() ) );
                    if( !ref_set.isEmpty() && ref_set.containsValue( Boolean.FALSE ) )
                    {
                        fileResult.add( ValidationMessage.error( "Unable to find reference sequence(s) from the CRAM reference registry: "
                                + ref_set.entrySet()
                                .stream()
                                .filter( e -> !e.getValue() )
                                .map( e -> e.getKey() )
                                .collect(Collectors.toList() ) ) );
                    }
                } catch( IOException ex )
                {
                    fileResult.add( ValidationMessage.error( ex ) );
                }
            }
        }

        readSamFile( fileResult, rf, paired, true );
    }
    
    
    public void readBamFile( ValidationResult result, RawReadsFile rf, AtomicBoolean paired ) throws IOException {
        readSamFile(result, rf, paired, false);
    }

    private void readSamFile( ValidationResult result, RawReadsFile rf, AtomicBoolean paired, boolean useCramRefSource )
            throws IOException {
        long read_no = 0;
        long reads_cnt = 0;

        try
        {
            result.add(ValidationMessage.info("Processing file"));

            File file = new File( rf.getFilename() );

            Log.setGlobalLogLevel( LogLevel.ERROR );

            SamReaderFactory.setDefaultValidationStringency( ValidationStringency.SILENT );
            SamReaderFactory factory = SamReaderFactory.make();
            factory.enable( SamReaderFactory.Option.DONT_MEMORY_MAP_INDEX );
            factory.validationStringency( ValidationStringency.SILENT );
            factory.samRecordFactory( DefaultSAMRecordFactory.getInstance() );

            if (useCramRefSource) {
                ENAReferenceSource reference_source = new ENAReferenceSource();
                reference_source.setLoggerWrapper( new ENAReferenceSource.LoggerWrapper()
                {
                    @Override public void
                    error( Object... messageParts )
                    {
                        result.add( ValidationMessage.error( null == messageParts ? "null" : String.valueOf( Arrays.asList( messageParts ) ) ) );
                    }

                    @Override public void
                    warn( Object... messageParts )
                    {
                        result.add( ValidationMessage.info( null == messageParts ? "null" : String.valueOf( Arrays.asList( messageParts ) ) ) );
                    }

                    @Override public void
                    info( Object... messageParts )
                    {
                        result.add( ValidationMessage.info( null == messageParts ? "null" : String.valueOf( Arrays.asList( messageParts ) ) ) );
                    }
                } );

                result.add( ValidationMessage.info( "REF_PATH  " + reference_source.getRefPathList() ) );
                result.add( ValidationMessage.info( "REF_CACHE " + reference_source.getRefCacheList() ) );

                factory.referenceSource( reference_source );
            }

//                SAMTextHeaderCodec codec = new SAMTextHeaderCodec();
//                codec.setValidationStringency( ValidationStringency.SILENT );

            SamInputResource ir = SamInputResource.of( file );
            File indexMaybe = SamFiles.findIndex( file );
            result.add( ValidationMessage.info( "proposed index: " + indexMaybe ) );

            if( null != indexMaybe )
                ir.index( indexMaybe );

            try( SamReader reader = factory.open(ir) )
            {
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
                        ValidationResult readResult = result.create( new ValidationOrigin( "read number", read_no ) );
                        readResult.add( ValidationMessage.error( "Mismatch between length of read bases and qualities" ) );
                    }

                    paired.compareAndSet( false, record.getReadPairedFlag() );
                    reads_cnt++;
                    if( 0 == reads_cnt % print_freq )
                        logProcessedReadNumber( reads_cnt );
                }

                logProcessedReadNumber( reads_cnt );

            }

            result.add( ValidationMessage.info( "LibraryLayout: " + ( paired.get() ? "PAIRED" : "SINGLE" ) ) );

            if( 0 == reads_cnt )
                result.add( ValidationMessage.error( "File contains no valid reads" ) );

        } catch( SAMFormatException | CRAMException e )
        {
            result.add( ValidationMessage.error( e ) );
        }
    }
}
