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

import org.apache.commons.compress.compressors.bzip2.BZip2CompressorInputStream;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import uk.ac.ebi.ena.readtools.common.reads.QualityNormalizer;
import uk.ac.ebi.ena.readtools.common.reads.normalizers.htsjdk.StandardQualityNormalizer;
import uk.ac.ebi.ena.readtools.loader.common.consumer.Spot;
import uk.ac.ebi.ena.readtools.loader.common.consumer.DataConsumer;
import uk.ac.ebi.ena.readtools.loader.common.consumer.DataConsumerException;
import uk.ac.ebi.ena.readtools.loader.common.producer.DataProducerException;
import uk.ac.ebi.ena.readtools.loader.common.producer.DataSpotProducer;
import uk.ac.ebi.ena.readtools.loader.fastq.DataSpot;
import uk.ac.ebi.ena.readtools.loader.fastq.FastqIterativeConsumer;
import uk.ac.ebi.ena.readtools.loader.fastq.FastqIterativeConsumer.READ_TYPE;
import uk.ac.ebi.ena.readtools.loader.fastq.FastqSpot;
import uk.ac.ebi.ena.readtools.loader.fastq.PairedFastqConsumer;
import uk.ac.ebi.ena.readtools.utils.Utils;
import uk.ac.ebi.ena.webin.cli.validator.message.ValidationMessage;
import uk.ac.ebi.ena.webin.cli.validator.message.ValidationOrigin;
import uk.ac.ebi.ena.webin.cli.validator.message.ValidationResult;

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


    private QualityNormalizer getQualityNormalizer(RawReadsFile rf ) {
        return Utils.getQualityNormalizer(Utils.detectFastqQualityFormat(rf.getFilename(), null));
    }
    
    
    private DataProducerException
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
            
            DataSpotProducer dp = new DataSpotProducer( is, normalizer, "" );
            dp.setName( stream_name );
            
            dp.setConsumer(new DataConsumer<DataSpot, Spot>()
            {
                @Override
                public void cascadeErrors() throws DataConsumerException { }

                @Override public void
                consume(DataSpot spot )
                {
                	String readKey;
                	String readIndex;

                	try
                	{
                		readKey = PairedFastqConsumer.getReadKey( spot.name);
                		readIndex = PairedFastqConsumer.getReadIndex( spot.name);
                	} catch ( DataConsumerException dee )
                	{
                    	readKey  = spot.name;
                    	readIndex = stream_name;
                	}
                	
                    if( labels.size() < MAX_LABEL_SET_SIZE )
                        labels.add( readIndex );
                    
                    count.incrementAndGet();
                    pairing.add( readKey );
                    duplications.add( spot.name);
                    
                    if( 0 == count.get() % print_freq )
                        logProcessedReadNumber( count.get() );
                }

                @Override
                public void setConsumer(DataConsumer<Spot, ?> dataConsumer) {
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

            if( !dp.isOk() && !( dp.getStoredException() instanceof DataProducerException) && !( dp.getStoredException() instanceof InvocationTargetException ) )
                throw dp.getStoredException();
            
            Throwable t = dp.isOk() ? dp.getReadRecordCount() > 0 ? null : new DataProducerException( 0, "Empty file" )
                                    : dp.getStoredException();
            if( dp.isOk() && null == t )
                return null;
            
            t = null == t ? new DataProducerException( -1, "Unknown failure" ) : t;
            t = t instanceof InvocationTargetException ? t.getCause() : t;
            t = t instanceof DataProducerException ? t : new DataProducerException( t );
            
            DataProducerException result = (DataProducerException) t;
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
        DataProducerException t = read( rf, labelset, pairing, duplications, count );
                        
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

            FastqIterativeConsumer wrapper = new FastqIterativeConsumer();
            wrapper.setFiles( new File[] { new File( rf.getFilename() ) } );  
            wrapper.setNormalizers( new QualityNormalizer[] { new StandardQualityNormalizer() } );
            wrapper.setReadType( READ_TYPE.SINGLE );

            Iterator<String> read_name_iterator = new DelegateIterator<FastqSpot, String>( wrapper.iterator() ) {
                @Override public String convert( FastqSpot obj )
                {
                    return obj.forward.name;
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
