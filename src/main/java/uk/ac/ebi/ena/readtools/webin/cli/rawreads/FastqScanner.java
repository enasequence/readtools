/*
 * Copyright 2018-2019 EMBL - European Bioinformatics Institute
 * Licensed under the Apache License, Version 2.0 (the "License"); you may not use this
 * file except in compliance with the License. You may obtain a copy of the License at
 * http://www.apache.org/licenses/LICENSE-2.0
 * Unless required by applicable law or agreed to in writing, software distributed under the
 * License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
 * CONDITIONS OF ANY KIND, either express or implied. See the License for the
 * specific language governing permissions and limitations under the License.
 */
package uk.ac.ebi.ena.readtools.webin.cli.rawreads;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.lang.reflect.InvocationTargetException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
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
import java.util.zip.GZIPInputStream;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import org.apache.commons.compress.compressors.bzip2.BZip2CompressorInputStream;

import uk.ac.ebi.ena.readtools.loader.common.QualityNormalizer;
import uk.ac.ebi.ena.readtools.loader.common.eater.DataEaterException;
import uk.ac.ebi.ena.readtools.loader.common.eater.NullDataEater;
import uk.ac.ebi.ena.readtools.loader.common.feeder.AbstractDataFeeder;
import uk.ac.ebi.ena.readtools.loader.common.feeder.DataFeederException;
import uk.ac.ebi.ena.readtools.loader.fastq.DataSpot;
import uk.ac.ebi.ena.readtools.loader.fastq.DataSpot.DataSpotParams;
import uk.ac.ebi.ena.readtools.loader.fastq.IlluminaIterativeEater;
import uk.ac.ebi.ena.readtools.loader.fastq.IlluminaIterativeEater.READ_TYPE;
import uk.ac.ebi.ena.readtools.loader.fastq.IlluminaPairedDataEater;
import uk.ac.ebi.ena.readtools.loader.fastq.IlluminaSpot;
import uk.ac.ebi.ena.webin.cli.validator.message.ValidationMessage;
import uk.ac.ebi.ena.webin.cli.validator.message.ValidationOrigin;
import uk.ac.ebi.ena.webin.cli.validator.message.ValidationResult;


public abstract class 
FastqScanner
{
    private static final int MAX_LABEL_SET_SIZE = 10;
    private static final int PAIRING_THRESHOLD = 20;
    
    //TODO remove duplication
    private static final int    DEFAULT_PRINT_FREQ = 1_000;
    private static final String PRINT_FREQ_PROPERTY_NAME = "webincli.scanner.print.freq";
    private static final int    print_freq = Integer.valueOf( System.getProperty( PRINT_FREQ_PROPERTY_NAME, String.valueOf( DEFAULT_PRINT_FREQ ) ) );    
    
    private final int expected_size;
    private final Set<String> labelset = new HashSet<>();
    private final AtomicBoolean paired = new AtomicBoolean();

    private static final Logger log = LoggerFactory.getLogger( FastqScanner.class );
    
    abstract protected void logProcessedReadNumber( long count );
    abstract protected void logFlushMsg( String message );
    
    
    public
    FastqScanner( int expected_size )
    {
        this.expected_size = expected_size;
    }

    
    private InputStream
    openFileInputStream( Path path )
    {
        final int marksize = 256;
        BufferedInputStream is;
        try 
        {
            is = new BufferedInputStream( Files.newInputStream( path ) );
            is.mark( marksize );
            try
            {
                return new BufferedInputStream( new GZIPInputStream( is ) );
            } catch( IOException gzip )
            {
                is.reset();
                try
                {
                    is.mark( marksize );
                    return new BufferedInputStream( new BZip2CompressorInputStream( is ) );
                } catch( IOException bzip )
                {
                    is.reset();
                    return is;
                }
            }
        } catch( IOException ex )
        {
            throw new RawReadsException( ex, ex.getMessage() );
        }
    }
    
    
    private QualityNormalizer
    getQualityNormalizer( RawReadsFile rf )
    {
        QualityNormalizer qn = QualityNormalizer.NONE;
        
        if( null != rf.getQualityScoringSystem() )
        {
            switch( rf.getQualityScoringSystem() )
            {
            default:
                throw new RawReadsException( "Scoring system: " + rf.getQualityScoringSystem() );
   
            case phred:
                switch( rf.getAsciiOffset() )
                {
                default:
                    throw new RawReadsException( "ASCII offset: " + rf.getAsciiOffset() );
                    
                case FROM33:
                    qn = QualityNormalizer.X;
                    break;
                    
                case FROM64:
                    qn = QualityNormalizer.X_2;
                    break;
                }
                break;
                
            case log_odds:
                qn = QualityNormalizer.SOLEXA;
                break;
            }
        }
        return qn;
    }
    
    
    private DataFeederException 
    read( RawReadsFile rf,
          Set<String>  labels,
          BloomWrapper pairing,
          BloomWrapper duplications,
          AtomicLong count ) throws Throwable
    {
        try( InputStream is = openFileInputStream( Paths.get( rf.getFilename() ) ) )
        {
            String stream_name = rf.getFilename();
            final QualityNormalizer normalizer = getQualityNormalizer( rf );
            
            AbstractDataFeeder<DataSpot> df = new AbstractDataFeeder<DataSpot>( is, DataSpot.class ) 
            {
                final DataSpotParams params = DataSpot.defaultParams();
                
                @Override protected DataSpot 
                newFeedable()
                {
                    return new DataSpot( normalizer, "", params );
                }
            };
            
            df.setName( stream_name );
            
            df.setEater( new NullDataEater<DataSpot>() 
            {
                @Override public void
                eat( DataSpot spot ) 
                {
                	String name  = null;
                	String label = null;
                			
                	try
                	{
                		name = IlluminaPairedDataEater.getReadnamePart( spot.bname, IlluminaPairedDataEater.KEY );
                		label = IlluminaPairedDataEater.getReadnamePart( spot.bname, IlluminaPairedDataEater.INDEX );
                	} catch ( DataEaterException dee )
                	{
                    	name  = spot.bname;
                    	label = stream_name;
                	}
                	
                    if( labels.size() < MAX_LABEL_SET_SIZE )
                        labels.add( label );
                    
                    count.incrementAndGet();
                    pairing.add( name );
                    duplications.add( spot.bname );
                    
                    if( 0 == count.get() % print_freq )
                        logProcessedReadNumber( count.get() );
                }  
            } );


            log.info( "Processing file " + rf.getFilename() );
            df.start();
            df.join();
            logProcessedReadNumber( count.get() );
            logFlushMsg( String.format( ", result: %s\n", null == df.getStoredException() ? "OK" : String.valueOf( df.getStoredException() ) ) );

            if( !df.isOk() && !( df.getStoredException() instanceof DataFeederException ) && !( df.getStoredException() instanceof InvocationTargetException ) )
                throw df.getStoredException();
            
            Throwable t = df.isOk() ? df.getFieldFeedCount() > 0 ? null : new DataFeederException( 0, "Empty file" ) 
                                    : df.getStoredException();
            if( df.isOk() && null == t )
                return null;
            
            t = null == t ? new DataFeederException( -1, "Unknown failure" ) : t;
            t = t instanceof InvocationTargetException ? t.getCause() : t;
            t = t instanceof DataFeederException ? t : new DataFeederException( t );
            
            DataFeederException result = (DataFeederException) t;
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
        DataFeederException t = read( rf, labelset, pairing, duplications, count );
                        
        if( null != t )
        {
            ValidationResult dataFeederResult = fileResult.create(new ValidationOrigin("line number", t.getLineNo()));
            dataFeederResult.add(ValidationMessage.error( t.getMessage()));
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
        
        BloomWrapper duplications = new BloomWrapper( expected_size );
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
                
        if( 2 == labelset.size() )
        {
            paired.set( true );
            double pairing_level = (double)pairing.getPossibleDuplicateCount() / ( (double)( pairing.getAddsNumber() - pairing.getPossibleDuplicateCount() ) );
            pairing_level = 100 * ( 1 < pairing_level ? 1/ pairing_level : pairing_level );
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

            long index = 1;

            IlluminaIterativeEater wrapper = new IlluminaIterativeEater();
            wrapper.setFiles( new File[] { new File( rf.getFilename() ) } );  
            wrapper.setNormalizers( new QualityNormalizer[] { QualityNormalizer.SANGER } );
            wrapper.setReadType( READ_TYPE.SINGLE );

            Iterator<String> read_name_iterator = new DelegateIterator<IlluminaSpot, String>( wrapper.iterator() ) {
                @Override public String convert( IlluminaSpot obj )
                {
                    return obj.read_name[ IlluminaSpot.FORWARD ];
                }
            };
            
            while( read_name_iterator.hasNext() )
            {
                String read_name = read_name_iterator.next();
                if( duplications.getSuspected().contains( read_name ) )
                {
                    counts.put( read_name, counts.getOrDefault( read_name, 0 ) + 1 );
                    Set<String> dlist = results.getOrDefault( read_name, new LinkedHashSet<>() );
                    dlist.add( String.format( "%s, read %s", rf.getFilename(), index ) );
                    results.put( read_name, dlist );
                }
                index ++;
            }
        }
        
        return results.entrySet()
                      .stream()
                      .filter( e-> counts.get(e.getKey()) > 1 )
                      .limit( limit )
                      .collect( Collectors.toMap( e -> e.getKey(), e -> e.getValue(), ( v1, v2 ) -> v1, LinkedHashMap::new ) );
    }
}
