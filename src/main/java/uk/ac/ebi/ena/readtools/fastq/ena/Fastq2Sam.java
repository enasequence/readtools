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
package uk.ac.ebi.ena.readtools.fastq.ena;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParameterException;
import com.beust.jcommander.Parameters;
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.util.FastqQualityFormat;
import htsjdk.samtools.util.QualityEncodingDetector;
import uk.ac.ebi.ena.readtools.common.reads.QualityNormalizer;
import uk.ac.ebi.ena.readtools.common.reads.normalizers.htsjdk.IlluminaQualityNormalizer;
import uk.ac.ebi.ena.readtools.common.reads.normalizers.htsjdk.SolexaQualityNormalizer;
import uk.ac.ebi.ena.readtools.common.reads.normalizers.htsjdk.StandardQualityNormalizer;
import uk.ac.ebi.ena.readtools.loader.common.FileCompression;
import uk.ac.ebi.ena.readtools.loader.common.consumer.DataConsumer;
import uk.ac.ebi.ena.readtools.loader.common.producer.AbstractDataProducer;
import uk.ac.ebi.ena.readtools.loader.common.producer.DataProducerException;
import uk.ac.ebi.ena.readtools.loader.fastq.DataSpot;
import uk.ac.ebi.ena.readtools.loader.fastq.DataSpot.DataSpotParams;
import uk.ac.ebi.ena.readtools.loader.fastq.PairedFastqConsumer;
import uk.ac.ebi.ena.readtools.loader.fastq.SingleFastqConsumer;
import uk.ac.ebi.ena.readtools.loader.fastq.FastqSpot;
import uk.ac.ebi.ena.readtools.utils.Utils;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class Fastq2Sam {

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
            new Fastq2Sam().create( p );
            System.exit( Params.OK_CODE );
            
        } catch( Throwable e ) {
            e.printStackTrace();
            System.exit( Params.FAILURE );
        }
    }
    
    public void create( Params p ) throws IOException {
        DataConsumer<DataSpot, FastqSpot> dataSpotToFastqSpotConsumer = null;

        if( null == p.files || p.files.size() < 1 || p.files.size() > 2) {
            throw new IllegalArgumentException("Invalid number of input files : " + p.files.size());
        } else if( 1 == p.files.size() ) {
            //single
            dataSpotToFastqSpotConsumer = (DataConsumer<DataSpot, FastqSpot>)new SingleFastqConsumer();
        } else if( 2 == p.files.size() ) {
            //same file names;
            if( p.files.get( 0 ).equals( p.files.get( 1 ) ) ) {
                throw new IllegalArgumentException(
                        "Paired files cannot be same. File1 : " + p.files.get(0) + ", File2 : " + p.files.get(1));
            }

            dataSpotToFastqSpotConsumer = (DataConsumer<DataSpot, FastqSpot>) new PairedFastqConsumer(
                    new File( p.tmp_root ), p.spill_page_size );
        }
        
        if( p.verbose ) {
            System.out.println( "Files to process: " );
            for( String f_name: p.files )
                System.out.println( " " + f_name );
        }

        FastqQualityFormat qualityFormat = Utils.detectFastqQualityFormat(p.files.get(0),
                p.files.size() == 2 ? p.files.get(0) : null);

        Fastq2BamConsumer fastqSpotToBamConsumer = new Fastq2BamConsumer(
                determineQualityNormalizer(qualityFormat), p.sample_name, p.data_file, p.tmp_root);
        
        dataSpotToFastqSpotConsumer.setConsumer( fastqSpotToBamConsumer );
        
        ArrayList<AbstractDataProducer<?>> producers = new ArrayList<AbstractDataProducer<?>>();

        int attr = 1;
        for( String f_name: p.files ) {
            final String default_attr = Integer.toString( attr ++ );

            AbstractDataProducer<DataSpot> producer = new AbstractDataProducer<DataSpot>(
                    FileCompression.valueOf( p.compression ).open( f_name, p.use_tar )) {
                final DataSpotParams params = DataSpot.defaultParams();
                
                @Override
                protected DataSpot
                newProducible()
                {
                    return new DataSpot( null, default_attr, params );
                }
            };

            producer.setConsumer( dataSpotToFastqSpotConsumer );
            producer.setName( f_name );
            producers.add( producer );
            producer.start();
        }
        
        boolean again = false;
        do {
            for( AbstractDataProducer<?> producer : producers ) {
                if( producer.isAlive() && producer.isOk() ) {
                    try {
                        producer.join();
                    } catch( InterruptedException ie ) {
                        again = true;
                    }
                } else if( !producer.isOk() ) {
                    throw new DataProducerException( producer.getStoredException() );
                }
            }
        } while( again );
       
        for( AbstractDataProducer<?> producer : producers ) {
            if( !producer.isOk() ) {
                throw new DataProducerException( producer.getStoredException() );
            }
        }

        dataSpotToFastqSpotConsumer.cascadeErrors();
        fastqSpotToBamConsumer.unwind();
        System.out.println( "DONE" );
    }

    private QualityNormalizer determineQualityNormalizer(FastqQualityFormat qualityType) {
        switch (qualityType)  {
            case Standard:
                return new StandardQualityNormalizer();
            case Solexa:
                return new SolexaQualityNormalizer();
            case Illumina:
                return new IlluminaQualityNormalizer();
            default:
                throw new IllegalArgumentException("Unexpected fastq quality format provided : " + qualityType);
        }
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
            return String.format( "CommonParams:\nfiles: %s\ncompression: %s\ndata_file: %s",
                    files,
                    compression,
                    data_file);
        }
    }
}
