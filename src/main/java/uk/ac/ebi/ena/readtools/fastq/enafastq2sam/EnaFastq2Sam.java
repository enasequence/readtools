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
package uk.ac.ebi.ena.readtools.fastq.enafastq2sam;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParameterException;
import com.beust.jcommander.Parameters;
import uk.ac.ebi.ena.readtools.loader.common.FileCompression;
import uk.ac.ebi.ena.readtools.loader.common.QualityConverter;
import uk.ac.ebi.ena.readtools.loader.common.QualityNormalizer;
import uk.ac.ebi.ena.readtools.loader.common.eater.DataConsumer;
import uk.ac.ebi.ena.readtools.loader.common.eater.DataConsumerException;
import uk.ac.ebi.ena.readtools.loader.common.feeder.AbstractDataProducer;
import uk.ac.ebi.ena.readtools.loader.common.feeder.DataProducerException;
import uk.ac.ebi.ena.readtools.loader.fastq.DataSpot;
import uk.ac.ebi.ena.readtools.loader.fastq.DataSpot.DataSpotParams;
import uk.ac.ebi.ena.readtools.loader.fastq.IlluminaPairedDataConsumer;
import uk.ac.ebi.ena.readtools.loader.fastq.IlluminaSingleDataConsumer;
import uk.ac.ebi.ena.readtools.loader.fastq.IlluminaSpot;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class EnaFastq2Sam {

    public static void main( String[] args ) {
        final Params p = new Params();
        JCommander jc = new JCommander( p );
        try {
            jc.parse( args );
        } catch( ParameterException pe ) {
            jc.usage();
            System.exit( Params.PARAM_ERROR );
        }
        
        if( p.help ) {
            jc.usage();
            System.exit( Params.OK_CODE );
        }
        
        try {
            new EnaFastq2Sam().create( p );
            System.exit( Params.OK_CODE );
            
        } catch( Throwable e ) {
            e.printStackTrace();
            System.exit( Params.FAILURE );
        }
    }
    
    public void create( Params p ) throws SecurityException, FileNotFoundException, DataProducerException,
            NoSuchMethodException, IOException, DataConsumerException {

        if( null == p.files || p.files.size() < 1 )
            throw new IllegalStateException( "No input files" );
        
        if( p.verbose ) {
            System.out.println( "Files to process: " );
            for( String f_name: p.files )
                System.out.println( " " + f_name );
        }
        
        ILLUMINA_SETTINGS type = null;
        if( p.quality_type != null 
            && p.quality_type.trim().length() > 0 )
            type = ILLUMINA_SETTINGS.valueOf( p.quality_type );
        else
            type = ILLUMINA_SETTINGS.valueOf( QualityNormalizer.valueOf( p.quality_encoding ),
                                              QualityConverter.valueOf( p.quality_scoring ) );
        
        DataConsumer<DataSpot, IlluminaSpot> spotToIlluminaSpotConsumer = null;
        
        if( 1 == p.files.size() ) {
            //single
            spotToIlluminaSpotConsumer = (DataConsumer<DataSpot, IlluminaSpot>)new IlluminaSingleDataConsumer();

            if( p.verbose )
                System.out.println( "Single file/read expected" );

        } else if( 2 == p.files.size() ) {
            //same file names;
            if( p.files.get( 0 ).equals( p.files.get( 1 ) ) ) {
                p.files.remove( 1 );
                
                if( p.verbose )
                    System.out.println( "Single file/paired read expected" );
            } else {
                if( p.verbose )
                    System.out.println( "Paired file/read expected" );
            }

            spotToIlluminaSpotConsumer = (DataConsumer<DataSpot, IlluminaSpot>) new IlluminaPairedDataConsumer( new File( p.tmp_root ), p.spill_page_size );
            
        } else {
            throw new IllegalArgumentException("More than two input files not allowed.");
        }

        IlluminaBamConsumer illuminaSpotToBamConsumer = new IlluminaBamConsumer(type, p.sample_name, p.data_file);
        
        spotToIlluminaSpotConsumer.setConsumer( illuminaSpotToBamConsumer );
        
        ArrayList<AbstractDataProducer<?>> producers = new ArrayList<AbstractDataProducer<?>>();
        
        final QualityNormalizer normalizer = type.getNormalizer();
        int attr = 1;
        for( String f_name: p.files ) {
            final String default_attr = Integer.toString( attr ++ );

            AbstractDataProducer<?> producer = new AbstractDataProducer<DataSpot>(
                    FileCompression.valueOf( p.compression ).open( f_name, p.use_tar ), DataSpot.class ) {
                final DataSpotParams params = DataSpot.defaultParams();
                
                @Override
                protected DataSpot
                newProducible()
                {
                    return new DataSpot( normalizer, default_attr, params );
                }
            }.setConsumer( spotToIlluminaSpotConsumer );
            
            producer.setName( f_name );
            producers.add( producer );
            producer.start();
        }
        
        boolean again = false;
        do {
            for( AbstractDataProducer<?> feeder : producers ) {
                if( feeder.isAlive() && feeder.isOk() ) {
                    try {
                        feeder.join();
                    } catch( InterruptedException ie ) {
                        again = true;
                    }
                } else if( !feeder.isOk() ) {
                    throw new DataProducerException( feeder.getStoredException() );
                }
            }
        } while( again );
       
        for( AbstractDataProducer<?> producer : producers ) {
            if( !producer.isOk() ) {
                throw new DataProducerException( producer.getStoredException() );
            }
        }

        spotToIlluminaSpotConsumer.cascadeErrors();
        illuminaSpotToBamConsumer.unwind();
        System.out.println( "DONE" );
    }

    @Parameters(commandDescription = "FastQ to SAM conversion.")
    public static class Params {

        final public static int PARAM_ERROR = 2;
        final public static int OK_CODE     = 0;
        final public static int FAILURE     = 1;

        @Parameter( names = { "-h", "--help" } )
        public boolean help = false;

        @Parameter( names = { "-f", "--file" }, description = "files to be loaded, repeat option and parameter in case of more than one files (NB: order matters!)" )
        public List<String> files;

        @Parameter( names = { "-c", "--compression"}, required = false, description = "compression (applied for all input files), supported values: BZ2, GZIP or GZ, ZIP, BGZIP or BGZ and NONE" )
        public String compression = FileCompression.NONE.name();

        @Parameter( names = { "-o", "--output-data-file" }, description = "Output file" )
        public String data_file = "data.tmp";

        @Parameter( names = { "-q", "--quality-type" }, description = "types: SANGER, SOLEXA, ILLUMINA_1_3, ILLUMINA_1_5, if set - overrides values of -qs and -qe" )
        public String quality_type = null;

        @Parameter( names = { "-qe", "--quality-encoding" }, description = "X|X_2" )
        public String quality_encoding = QualityNormalizer.X.toString();

        @Parameter( names = { "-qs", "--quality-scoring" }, description = "PHRED|LOGODDS" )
        public String quality_scoring = QualityConverter.PHRED.toString();

        @Parameter( names = { "-v", "--verbose" }, description = "Verbose" )
        public boolean verbose = false;

        @Parameter( names = { "-tar" }, description = "Set this flag if input files are tarred" )
        public boolean use_tar = false;

        @Parameter( names = { "-sps", "-spill-page-size" }, description = "Spill page size, depends on maximum of avaliable memory and size of unassembled record pool" )
        public int spill_page_size = 7000000;

        @Parameter( names = { "-tmp", "--tmp-root" }, description = "Folder to store temporary output in case of paired reads assembly. Should be large enough" )
        public String tmp_root = ".";

        @Parameter( names = { "-sm", "--sample-name" }, required = true, description = "Value to use for SAM header SM. Required." )
        public String sample_name = null;

        public String toString() {
            return String.format( "CommonParams:\nfiles: %s\ncompression: %s\ndata_file: %s\nquality_type: %s\nquality_encoding: %s\nquality_scoring: %s",
                    files,
                    compression,
                    data_file,
                    quality_type,
                    quality_encoding,
                    quality_scoring );
        }
    }
}
