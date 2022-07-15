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
package uk.ac.ebi.ena.readtools.webin.cli.rawreads;

import java.io.File;
import java.io.InputStream;
import java.lang.reflect.InvocationTargetException;
import java.nio.file.Paths;
import java.time.Duration;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicLong;
import java.util.stream.Collectors;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import uk.ac.ebi.ena.readtools.common.reads.QualityNormalizer;
import uk.ac.ebi.ena.readtools.common.reads.normalizers.htsjdk.StandardQualityNormalizer;
import uk.ac.ebi.ena.readtools.loader.common.converter.ConverterException;
import uk.ac.ebi.ena.readtools.loader.common.converter.FastqReadReadConverter;
import uk.ac.ebi.ena.readtools.loader.common.writer.ReadWriter;
import uk.ac.ebi.ena.readtools.loader.common.writer.ReadWriterException;
import uk.ac.ebi.ena.readtools.loader.common.writer.Spot;
import uk.ac.ebi.ena.readtools.loader.fastq.FastqIterativeWriter;
import uk.ac.ebi.ena.readtools.loader.fastq.FastqIterativeWriter.READ_TYPE;
import uk.ac.ebi.ena.readtools.loader.fastq.PairedFastqWriter;
import uk.ac.ebi.ena.readtools.loader.fastq.PairedRead;
import uk.ac.ebi.ena.readtools.loader.fastq.Read;
import uk.ac.ebi.ena.readtools.utils.Utils;
import uk.ac.ebi.ena.webin.cli.validator.message.ValidationMessage;
import uk.ac.ebi.ena.webin.cli.validator.message.ValidationOrigin;
import uk.ac.ebi.ena.webin.cli.validator.message.ValidationResult;


public abstract class 
FastqScanner
{
    private static final int MAX_LABEL_SET_SIZE = 10;

    /**
     * If the files are paired perfectly and there are no invalid duplicates then pairing level will be above this threshold
     * and equal to 100%.<br><br>
     * If there far are too many duplicates present then pairing level will result in more than 100%, the pairing calculation
     * algorithm will identify this suspicious increase and lower the pairing percentage in such a way that if there are
     * extremely large number of duplicates then pairing percentage will fall below this threshold. This then will be a sign
     * that files are either malformed (have lots of invalid duplicates) or not properly paired.<br><br>
     * To summarise, the calculated pairing percentage will be lower than this threshold when either:<br>
     * 1) There is a large number of duplicate reads in files<br>
     * 2) The reads files are incorrectly paired i.e. large number of unpaired reads than paired reads.
     */
    private static final int PAIRING_THRESHOLD = 20;
    
    //TODO remove duplication
    private static final int    DEFAULT_PRINT_FREQ = 1_000;
    private static final String PRINT_FREQ_PROPERTY_NAME = "webincli.scanner.print.freq";
    private static final int    print_freq = Integer.valueOf( System.getProperty( PRINT_FREQ_PROPERTY_NAME, String.valueOf( DEFAULT_PRINT_FREQ ) ) );

    private static final Logger log = LoggerFactory.getLogger( FastqScanner.class );
    
    private final int expected_size;
    private final Set<String> labelset = new HashSet<>();
    private final AtomicBoolean paired = new AtomicBoolean();

    private final Duration runDuration;
    
    abstract protected void logProcessedReadNumber( long count );
    abstract protected void logFlushMsg( String message );
    
    public
    FastqScanner( int expected_size )
    {
        this(null, expected_size);
    }

    public FastqScanner(Duration runDuration, int expected_size )
    {
        this.runDuration = runDuration;
        this.expected_size = expected_size;
    }
    
    
    private ConverterException
    read( RawReadsFile rf,
          Set<String>  labels,
          BloomWrapper pairing,
          BloomWrapper duplications,
          AtomicLong count ) throws Throwable
    {
        try( InputStream is = Utils.openFastqInputStream( Paths.get( rf.getFilename() ) ) )
        {
            String stream_name = rf.getFilename();

            FastqReadReadConverter dp = new FastqReadReadConverter(
                    is, runDuration, "", rf.getFilename());
            dp.setName( stream_name );
            
            dp.setWriter(new ReadWriter<Read, Spot>()
            {
                @Override
                public void cascadeErrors() throws ReadWriterException { }

                @Override public void
                write(Read spot )
                {
                	String readNameWithoutIndex;
                	String readIndex;

                	try
                	{
                		readNameWithoutIndex = PairedFastqWriter.getReadKey( spot.name);
                		readIndex = PairedFastqWriter.getReadIndex( spot.name);
                	} catch ( ReadWriterException dee )
                	{
                    	readNameWithoutIndex  = spot.name;
                    	readIndex = stream_name;
                	}
                	
                    if( labels.size() < MAX_LABEL_SET_SIZE ) {
                        labels.add(readIndex);
                    }
                    
                    count.incrementAndGet();
                    pairing.add( readNameWithoutIndex );
                    duplications.add( spot.name);
                    
                    if( 0 == count.get() % print_freq )
                        logProcessedReadNumber( count.get() );
                }

                @Override
                public void setWriter(ReadWriter<Spot, ?> readWriter) {
                    throw new RuntimeException( "Not implemented" );
                }

                @Override
                public boolean isOk() {
                    return true;
                }
            } );


            log.info( "Processing file " + rf.getFilename() );
            dp.start();
            dp.join();
            logProcessedReadNumber( count.get() );
            logFlushMsg( String.format( ", result: %s\n", null == dp.getStoredException() ? "OK" : String.valueOf( dp.getStoredException() ) ) );

            if( !dp.isOk() && !( dp.getStoredException() instanceof ConverterException) && !( dp.getStoredException() instanceof InvocationTargetException ) )
                throw dp.getStoredException();
            
            Throwable t = dp.isOk() ? dp.getReadCount() > 0 ? null : new ConverterException( 0, "Empty file" )
                                    : dp.getStoredException();
            if( dp.isOk() && null == t )
                return null;
            
            t = null == t ? new ConverterException( -1, "Unknown failure" ) : t;
            t = t instanceof InvocationTargetException ? t.getCause() : t;
            t = t instanceof ConverterException ? t : new ConverterException( t );
            
            ConverterException result = (ConverterException) t;
            return result;
        }
    }

    
    private void
    checkSingleFile( ValidationResult fileResult,
                     RawReadsFile rf,
                     Set<String> labelset,
                     BloomWrapper pairing,
                     BloomWrapper duplications ) throws Throwable
    {
        AtomicLong count = new AtomicLong();
        ConverterException converterException = read( rf, labelset, pairing, duplications, count );
                        
        if( null != converterException)
        {
            ValidationMessage dataProducerError = ValidationMessage.error( converterException.getMessage() );
            dataProducerError.appendOrigin(new ValidationOrigin("line number", converterException.getLineNo()));
            fileResult.add(dataProducerError);
        } else
        {
            fileResult.add(ValidationMessage.info( String.format( "Collected %d reads", count.get())));
            fileResult.add(ValidationMessage.info( String.format( "Collected %d read labels: %s", labelset.size(), labelset )));
            fileResult.add(ValidationMessage.info( String.format( "Has possible duplicate(s): " + duplications.hasPossibleDuplicates())));
        }
    }
    
    public boolean 
    getPaired()
    {
        return this.paired.get();
    }
    
    public void
    checkFiles( ValidationResult result, RawReadsFile... rfs ) throws Throwable
    {
        if( null == rfs || rfs.length != 1 && rfs.length != 2 )
        {
            //terminal error
            result.add( ValidationMessage.error( "Unusual amount of files. Can accept only 1 or 2, but got " + ( null == rfs ? "null" : rfs.length ) ) );
        }

        /** Should ideally have low to 0 number of duplicates. */
        BloomWrapper duplications = new BloomWrapper( expected_size );

        /** Should ideally have high but not too high number of duplicates. A very high number will point to presence
         * of invalid duplicates.
         * To keep memory consumption lower, because we can tolerate false positive here, use lower expected read size. */
        BloomWrapper pairing = new BloomWrapper( expected_size / 10 );

        for( RawReadsFile rf : rfs )
        {
            ValidationResult fileResult;
            if( rf.getReportFile() == null )
            {
                fileResult = result.create(new ValidationOrigin("file", rf.getFilename()));
            } else
            {
                fileResult = result.create( rf.getReportFile().toFile(), new ValidationOrigin("file", rf.getFilename() ) );
            }

            Set<String> flabelset = new HashSet<>();
            checkSingleFile( fileResult, rf, flabelset, pairing, duplications );
            labelset.addAll( flabelset );
        }
        
        if (!result.isValid())
            return;

        //following is only relevant if we are dealing with 2 fastqs which are supposedly paired.
        if (rfs.length == 2) {
            if( 2 == labelset.size() )
            {
                paired.set( true );
                double pairing_level = calculatePairingLevel(pairing);
                result.add(ValidationMessage.info(String.format( "Pairing percentage: %.2f%%", pairing_level )));

                //TODO: estimate bloom false positives impact
                if( (double)PAIRING_THRESHOLD > pairing_level )
                {
                    //terminal error
                    result.add(ValidationMessage.error(String.format( "Detected paired fastq submission with less than %d%% of paired reads", PAIRING_THRESHOLD )));
                }
            } else if( labelset.size() > 2 )
            {
                result.add(ValidationMessage.error(String.format(
                        "When submitting paired reads using two Fastq files the reads must follow Illumina paired read naming conventions. "
                                + "This was not the case for the submitted Fastq files: %s. Unable to determine pairing from set: %s",
                        rfs,
                        labelset.stream().limit( 10 ).collect( Collectors.joining( ",", "", 10 < labelset.size() ? "..." : "" ) ) )));
            }
        }
        
        //extra check for suspected reads
        if( duplications.hasPossibleDuplicates() )
        {
            ValidationResult duplicationResult = result.create();

            // read name, list
            Map<String, Set<String>> duplicates = findAllduplications( duplications, 100, rfs );

            duplicates.entrySet().stream().forEach( e -> duplicationResult.add(ValidationMessage.error(
                String.format( "Multiple (%d) occurrences of read name \"%s\" at: %s\n",
                         e.getValue().size(),
                         e.getKey(),
                         e.getValue().toString()))));

            if( duplicationResult.isValid() ) {
                result.add(ValidationMessage.info( "No duplicate read names found." ));
            }
        }
    }

    
    private Map<String, Set<String>>
    findAllduplications(BloomWrapper duplications, int limit, RawReadsFile... rfs)
    {
        Map<String, Integer> counts = new HashMap<>( limit );
        Map<String, Set<String>> results = new LinkedHashMap<>( limit );
        for( RawReadsFile rf: rfs )
        {
            String msg = "Performing additional checks for file " + rf.getFilename();
            log.info( msg );

            FastqIterativeWriter wrapper = new FastqIterativeWriter();
            wrapper.setFiles( new File[] { new File( rf.getFilename() ) } );
            wrapper.setReadType( READ_TYPE.SINGLE );

            Iterator<String> read_name_iterator = new DelegateIterator<PairedRead, String>( wrapper.iterator() ) {
                @Override public String convert( PairedRead obj )
                {
                    return obj.forward.name;
                }
            };

            long index = 1;

            while( read_name_iterator.hasNext() )
            {
                String read_name = read_name_iterator.next();
                if( duplications.getPossibleDuplicates().contains( read_name ) )
                {
                    counts.put( read_name, counts.getOrDefault( read_name, 0 ) + 1 );
                    Set<String> dlist = results.getOrDefault( read_name, new LinkedHashSet<>() );
                    dlist.add( rf.getFilename() + ", read " + index);
                    results.put( read_name, dlist );
                }
                index ++;
            }
        }
        
        return results.entrySet()
                      .stream()
                      //only read names occurring more than once are considered duplicates
                      .filter( e-> counts.get(e.getKey()) > 1 )
                      .limit( limit )
                      .collect( Collectors.toMap( e -> e.getKey(), e -> e.getValue(), ( v1, v2 ) -> v1, LinkedHashMap::new ) );
    }

    /**
     * Best explained with a help of an example where AddCount = 10.<br><br>
     * AddCount : 10, PossibleDuplicates : 0, PairingLevel : 0.00%<br>
     * AddCount : 10, PossibleDuplicates : 1, PairingLevel : 11.11%<br>
     * AddCount : 10, PossibleDuplicates : 2, PairingLevel : 25.00%<br>
     * AddCount : 10, PossibleDuplicates : 3, PairingLevel : 42.86%<br>
     * AddCount : 10, PossibleDuplicates : 4, PairingLevel : 66.67%<br>
     * AddCount : 10, PossibleDuplicates : 5, PairingLevel : 100.00%<br>
     * (The further down we get from this point means the higher gets the number of invalid duplicate reads.)<br>
     * AddCount : 10, PossibleDuplicates : 6, PairingLevel : 66.67%<br>
     * AddCount : 10, PossibleDuplicates : 7, PairingLevel : 42.86%<br>
     * AddCount : 10, PossibleDuplicates : 8, PairingLevel : 25.00%<br>
     * AddCount : 10, PossibleDuplicates : 9, PairingLevel : 11.11%<br>
     * AddCount : 10, PossibleDuplicates : 10, PairingLevel : 0.00%<br>
     */
    private double calculatePairingLevel(BloomWrapper pairing) {
        double possiblePairedCount = pairing.getPossibleDuplicateCount();

        double addCount = pairing.getAddCount();

        double unpairedCount = addCount - possiblePairedCount;

        double pairedToUnpairedRatio = possiblePairedCount / unpairedCount;

        double pairingLevel;

        //Ratio being more than 1 indicates there are might be more than 1 duplicates (i.e. invalid duplicates) present
        //for some reads.
        if (pairedToUnpairedRatio > 1) {
            //Lower the pairing level in such case to so that if there are too many invalid duplicates the pairing level will
            //then go below the PAIRING_THRESHOLD and the validation will fail.
            pairingLevel = 1 / pairedToUnpairedRatio;
        } else {
            pairingLevel = pairedToUnpairedRatio;
        }

        pairingLevel = 100 * pairingLevel;

        return pairingLevel;
    }
}
