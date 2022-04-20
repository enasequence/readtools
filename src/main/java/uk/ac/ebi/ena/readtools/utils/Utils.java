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
package uk.ac.ebi.ena.readtools.utils;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;

import org.apache.commons.compress.compressors.bzip2.BZip2CompressorInputStream;

import htsjdk.samtools.fastq.AsyncFastqWriter;
import htsjdk.samtools.fastq.BasicFastqWriter;
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqWriter;
import htsjdk.samtools.util.FastqQualityFormat;
import htsjdk.samtools.util.QualityEncodingDetector;

import uk.ac.ebi.ena.readtools.common.reads.QualityNormalizer;
import uk.ac.ebi.ena.readtools.common.reads.normalizers.htsjdk.IlluminaQualityNormalizer;
import uk.ac.ebi.ena.readtools.common.reads.normalizers.htsjdk.SolexaQualityNormalizer;
import uk.ac.ebi.ena.readtools.common.reads.normalizers.htsjdk.StandardQualityNormalizer;
import uk.ac.ebi.ena.readtools.webin.cli.rawreads.RawReadsException;

public class Utils {

    private static final Pattern URACIL_PATTERN = Pattern.compile("u|U");

    /**
     * Generates a new Fastq file that has all of the U bases inside it replaced with T ones.
     *
     * @param inputFastq - Path to the input fastq file.
     * @param outputFastq - Path to the output fastq file. File extension determines compression.
     * @throws IOException
     */
    public static void replaceUracilBasesInFastq(String inputFastq, String outputFastq) throws IOException {
        Matcher matcher = URACIL_PATTERN.matcher("");

        File inp = new File(inputFastq);
        File out = new File(outputFastq);

        FastqWriter writer = new AsyncFastqWriter(new BasicFastqWriter(out), AsyncFastqWriter.DEFAULT_QUEUE_SIZE);

        FastqReader reader = new FastqReader(inp);
        for (FastqRecord record : reader) {
            String result = record.getReadString();

            if (matcher.reset(result).find()) {
                String replacement = matcher.group().equals("U") ? "T" : "t";

                result = matcher.replaceAll(replacement);
            }

            writer.write(new FastqRecord(record.getReadName(), result, record.getBaseQualityHeader(), record.getBaseQualityString()));
        }

        reader.close();
        writer.close();
    }

    /**
     * Replace U bases with T.
     *
     * @param bases
     * @return a new string of bases in which U bases are replaced with T's. Returns the input string as is otherwise.
     */
    public static String replaceUracilBases(String bases) {
        Matcher matcher = URACIL_PATTERN.matcher(bases);
        if (matcher.find()) {
            String replacement = matcher.group().equals("U") ? "T" : "t";

            return matcher.replaceAll(replacement);
        }

        return bases;
    }

    /**
     * Detects Fastq quality format by examining the given files. Second argument is optional but can be useful
     * when dealing with paired files.
     *
     * @param fastqFile1 - Path to fastq file. Cannot be null.
     * @param fastqFile2 - Path to fastq file. Optional.
     * @return
     */
    public static FastqQualityFormat detectFastqQualityFormat(String fastqFile1, String fastqFile2) {
        FastqReader reader1 = new FastqReader(null, new BufferedReader(new InputStreamReader(
                openFastqInputStream(Paths.get(fastqFile1)), StandardCharsets.UTF_8)), true);

        FastqReader reader2 = fastqFile2 == null ? null : new FastqReader(null, new BufferedReader(new InputStreamReader(
                openFastqInputStream(Paths.get(fastqFile2)), StandardCharsets.UTF_8)), true);

        final QualityEncodingDetector detector = new QualityEncodingDetector();

        if (reader2 == null) {
            detector.add(QualityEncodingDetector.DEFAULT_MAX_RECORDS_TO_ITERATE, reader1);
        } else {
            detector.add(QualityEncodingDetector.DEFAULT_MAX_RECORDS_TO_ITERATE, reader1, reader2);
            reader2.close();
        }

        reader1.close();

        final FastqQualityFormat qualityFormat =  detector.generateBestGuess(
                QualityEncodingDetector.FileContext.FASTQ, null);

        return qualityFormat;
    }

    public static QualityNormalizer getQualityNormalizer(FastqQualityFormat qualityType) {
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

    /**
     * Opens an input stream to Fastq file. Supports reading from *.gz and *.bzip2 archive formats.
     *
     * @param path
     * @return
     */
    public static InputStream openFastqInputStream(Path path ) {
        final int marksize = 256;
        BufferedInputStream is;

        try {
            is = new BufferedInputStream( Files.newInputStream( path ) );
            is.mark( marksize );

            try {
                return new BufferedInputStream( new GZIPInputStream( is ) );
            } catch( IOException gzip ) {
                is.reset();
                try {
                    is.mark( marksize );
                    return new BufferedInputStream( new BZip2CompressorInputStream( is, true ) );
                } catch( IOException bzip ) {
                    is.reset();
                    return is;
                }
            }
        } catch( IOException ex ) {
            throw new RawReadsException( ex, ex.getMessage() );
        }
    }
}
