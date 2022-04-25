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
package uk.ac.ebi.ena.readtools.fastq.ena;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParameterException;
import com.beust.jcommander.Parameters;
import htsjdk.samtools.util.FastqQualityFormat;
import uk.ac.ebi.ena.readtools.common.reads.QualityNormalizer;
import uk.ac.ebi.ena.readtools.loader.common.FileCompression;
import uk.ac.ebi.ena.readtools.loader.common.converter.ConverterException;
import uk.ac.ebi.ena.readtools.loader.common.converter.FastqReadConverter;
import uk.ac.ebi.ena.readtools.loader.common.writer.ReadWriter;
import uk.ac.ebi.ena.readtools.loader.common.writer.ReadWriterMemoryLimitException;
import uk.ac.ebi.ena.readtools.loader.fastq.PairedFastqWriter;
import uk.ac.ebi.ena.readtools.loader.fastq.PairedRead;
import uk.ac.ebi.ena.readtools.loader.fastq.Read;
import uk.ac.ebi.ena.readtools.loader.fastq.SingleFastqWriter;
import uk.ac.ebi.ena.readtools.utils.Utils;

import java.io.File;
import java.io.IOException;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

public class Sra2Sam {
    private long totalReadCount = 0, totalBaseCount = 0;

    public static void main(String[] args) {
        final Params p = new Params();
        JCommander jc = new JCommander(p);
        try {
            jc.parse(args);
        } catch (ParameterException pe) {
            jc.usage();
            System.exit(Params.PARAM_ERROR);
        }

        if (p.help) {
            jc.usage();
            System.exit(Params.OK_CODE);
        }

        try {
            new Sra2Sam().create(p);
            System.exit(Params.OK_CODE);

        } catch (Throwable e) {
            e.printStackTrace();
            System.exit(Params.FAILURE);
        }
    }

    public void create(Params p) throws IOException {
        ReadWriter<Read, PairedRead> dataSpotToFastqReadWriter = null;
        boolean paired = false;

        if (null == p.files || p.files.size() < 1 || p.files.size() > 2) {
            throw new IllegalArgumentException("Invalid number of input files : " + p.files.size());
        } else if (1 == p.files.size()) {
            //single
            dataSpotToFastqReadWriter = new SingleFastqWriter();
            paired = false;
        } else if (2 == p.files.size()) {
            //same file names;
            if (p.files.get(0).equals(p.files.get(1))) {
                throw new IllegalArgumentException(
                        "Paired files cannot be same. File1 : " + p.files.get(0) + ", File2 : " + p.files.get(1));
            }

            dataSpotToFastqReadWriter = new PairedFastqWriter(
                    new File(p.tmp_root), p.spill_page_size, p.spill_page_size_bytes, p.spill_abandon_limit_bytes);
            paired = true;
        }

        if (p.verbose) {
            System.out.println("Files to process: ");
            for (String f_name : p.files)
                System.out.println(" " + f_name);
        }

        FastqQualityFormat qualityFormat = Utils.detectFastqQualityFormat(
                p.files.get(0),
                p.files.size() == 2 ? p.files.get(1) : null);

        QualityNormalizer normalizer = Utils.getQualityNormalizer(qualityFormat);

        Sra2BamWriter fastqSpotToBamWriter = new Sra2BamWriter(
                normalizer, p.sample_name, p.data_file, p.tmp_root, p.convertUracil, paired);

        dataSpotToFastqReadWriter.setWriter(fastqSpotToBamWriter);

        ArrayList<FastqReadConverter> converters = new ArrayList<>();

        ExecutorService executorService = Executors.newCachedThreadPool();
        List<Future<?>> futures = Collections.synchronizedList(new ArrayList<>());

        int attr = 1;
        for (String f_name : p.files) {
            final String default_attr = Integer.toString(attr++);

            FastqReadConverter converter = new FastqReadConverter(
                    new Sra2SamInputStream(Paths.get(f_name)), default_attr, f_name);
            converter.setWriter(dataSpotToFastqReadWriter);
            converters.add(converter);

            futures.add(executorService.submit(converter::run));
        }

        futures.forEach(f -> {
            try {
                f.get();
            } catch (InterruptedException ignored) {
            } catch (ExecutionException e) {
                throw new ConverterException(e);
            }
        });

        for (FastqReadConverter converter : converters) {
            if (!converter.isOk()) {
                if (converter.getStoredException() instanceof ReadWriterMemoryLimitException) {
                    throw new ReadWriterMemoryLimitException(converter.getStoredException());
                }
                throw new RuntimeException(converter.getStoredException());
            }

            totalReadCount += converter.getReadCount();
            totalBaseCount += converter.getBaseCount();
        }

        dataSpotToFastqReadWriter.cascadeErrors();
        fastqSpotToBamWriter.unwind();
        System.out.println(String.format("READS: %d; BASES: %d", totalReadCount, totalBaseCount));
    }

    /**
     * @return Total number of reads that were processed.
     */
    public long getTotalReadCount() {
        return totalReadCount;
    }

    /**
     * @return Total number of bases accumulated across all reads.
     */
    public long getTotalBaseCount() {
        return totalBaseCount;
    }

    @Parameters(commandDescription = "FastQ to SAM conversion.")
    public static class Params {

        final public static int PARAM_ERROR = 2;
        final public static int OK_CODE = 0;
        final public static int FAILURE = 1;

        @Parameter(names = {"-h", "--help"})
        public boolean help = false;

        @Parameter(names = {"-f", "--file"}, description = "files to be loaded, repeat option and parameter in case of more than one files (NB: order matters!)")
        public List<String> files;

        @Parameter(names = {"-c", "--compression"}, required = false, description = "compression (applied for all input files), supported values: BZ2, GZIP or GZ, ZIP, BGZIP or BGZ and NONE")
        public String compression = FileCompression.NONE.name();

        @Parameter(names = {"-o", "--output-data-file"}, description = "Output file")
        public String data_file = "data.tmp";

        @Parameter(names = {"-v", "--verbose"}, description = "Verbose")
        public boolean verbose = false;

        @Parameter(names = {"-tar"}, description = "Set this flag if input files are tarred")
        public boolean use_tar = false;

        @Parameter(names = {"-sps", "-spill-page-size"}, description = "Spill page size, depends on maximum of available memory and size of un-assembled record pool")
        public int spill_page_size = 4_500_000;

        @Parameter(names = {"-sps", "-spill-page-size-bytes"}, description = "Spill page size in bytes, depends on maximum of available memory and size of un-assembled record pool")
        public long spill_page_size_bytes = 4L * 1024L * 1024L * 1024L;

        @Parameter(names = {"-sps", "-spill_abandon_limit_bytes"}, description = "Spill memory limit in bytes, processing fails when temp files total size reaches this limit, 0 == no limit")
        public long spill_abandon_limit_bytes = 10L * 1024L * 1024L * 1024L;

        @Parameter(names = {"-tmp", "--tmp-root"}, description = "Folder to store temporary output in case of paired reads assembly. Should be large enough")
        public String tmp_root = ".";

        @Parameter(names = {"-sm", "--sample-name"}, required = true, description = "Value to use for SAM header SM. Required.")
        public String sample_name = null;

        @Parameter(names = {"--convert-uracil"}, description = "Whether or not to convert Uracil bases [U, u] to [T, t]. Default is false.")
        public boolean convertUracil = false;

        public String toString() {
            return String.format("CommonParams:\nfiles: %s\ncompression: %s\ndata_file: %s",
                    files,
                    compression,
                    data_file);
        }
    }
}
