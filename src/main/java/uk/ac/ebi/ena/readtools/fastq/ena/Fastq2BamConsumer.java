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

import htsjdk.samtools.ReservedTagConstants;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import uk.ac.ebi.ena.readtools.common.reads.QualityNormalizer;
import uk.ac.ebi.ena.readtools.loader.common.InvalidBaseCharacterException;
import uk.ac.ebi.ena.readtools.loader.common.consumer.Spot;
import uk.ac.ebi.ena.readtools.loader.common.consumer.DataConsumer;
import uk.ac.ebi.ena.readtools.loader.common.consumer.DataConsumerException;
import uk.ac.ebi.ena.readtools.loader.fastq.DataSpot;
import uk.ac.ebi.ena.readtools.loader.fastq.FastqSpot;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Paths;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

/**
 * Accepts Fastq spot data and writes them out to a BAM file.
 */
public class Fastq2BamConsumer implements DataConsumer<FastqSpot, Spot> {

    private static final String DEFAULT_READ_GROUP_NAME = "A";

    private static final String VALID_DNA_CHARSET = "^[.acmgrsvtwyhkdbnNACMGRSVTWYHKDBN]+$";
    private static final Pattern VALID_DNA_CHARSET_PATTERN = Pattern.compile(VALID_DNA_CHARSET);

    private final QualityNormalizer qualityNormalizer;

    private final String sampleName;

    private final SAMFileWriter writer;

    private volatile boolean isOk = true;

    public Fastq2BamConsumer(QualityNormalizer qualityNormalizer, String sampleName, String outputFilePath, String tempDir) {
        this.qualityNormalizer = qualityNormalizer;
        this.sampleName = sampleName;

        if (sampleName == null || sampleName.trim().isEmpty()) {
            throw new IllegalArgumentException("Sample name is either null or empty.");
        }

        writer = new SAMFileWriterFactory().setTempDirectory(new File(tempDir)).makeSAMOrBAMWriter(
                createHeader(), false, Paths.get(outputFilePath));
    }

    @Override
    public void consume(FastqSpot spot) throws DataConsumerException {
        try {
            validate(spot);

            if (spot.isPaired()) {
                SAMRecord rec1 = createSamRecord(
                        true, spot.name, spot.forward.bases, spot.forward.quals);
                rec1.setFirstOfPairFlag(true);
                rec1.setSecondOfPairFlag(false);
                writer.addAlignment(rec1);

                SAMRecord rec2 = createSamRecord(
                        true, spot.name, spot.reverse.bases, spot.reverse.quals);
                rec2.setFirstOfPairFlag(false);
                rec2.setSecondOfPairFlag(true);
                writer.addAlignment(rec2);

            } else {
                DataSpot unpaired = spot.getUnpaired();
                SAMRecord rec = createSamRecord(false, spot.name, unpaired.bases, unpaired.quals);
                rec.setReadPairedFlag(false);
                writer.addAlignment(rec);
            }
        } catch (InvalidBaseCharacterException ex) {
            isOk = false;
            throw ex;
        } catch (Exception ex) {
            isOk = false;
            throw new DataConsumerException(ex);
        }
    }

    @Override
    public void cascadeErrors() throws DataConsumerException { }

    @Override
    public void setConsumer(DataConsumer<Spot, ? extends Spot> dataConsumer) { }

    @Override
    public boolean isOk() {
        return isOk;
    }

    public void unwind() {
        writer.close();
    }

    private SAMFileHeader createHeader() {
        final SAMReadGroupRecord rgroup = new SAMReadGroupRecord(DEFAULT_READ_GROUP_NAME);
        rgroup.setSample(this.sampleName);

        final SAMFileHeader header = new SAMFileHeader();
        header.addReadGroup(rgroup);
        header.setSortOrder(SAMFileHeader.SortOrder.queryname);

        return header;
    }

    private void validate(FastqSpot fastqSpot) {
        if (fastqSpot.forward != null) {
            Matcher matcher = VALID_DNA_CHARSET_PATTERN.matcher(fastqSpot.forward.bases);
            if( !matcher.matches() ) {
                handleInvalidDnaCharset(fastqSpot.forward.bases, matcher);
            }

            if( fastqSpot.forward.bases.length() != fastqSpot.forward.quals.length() )
                throw new IllegalArgumentException( String.format( "FATAL: Spot bases and qualities length do not match. Malformed spot\n%s\n", fastqSpot ) );
        }

        if (fastqSpot.reverse != null) {
            Matcher matcher = VALID_DNA_CHARSET_PATTERN.matcher(fastqSpot.reverse.bases);
            if( !matcher.matches() )
                handleInvalidDnaCharset(fastqSpot.reverse.bases, matcher);

            if( fastqSpot.reverse.bases.length() != fastqSpot.reverse.quals.length() )
                throw new IllegalArgumentException( String.format( "FATAL: Spot bases and qualities length do not match. Malformed spot\n%s\n", fastqSpot ) );
        }
    }

    private void handleInvalidDnaCharset(String bases, Matcher matcher) {
        Set<Character> invalidBaseChars = bases.chars()
                .mapToObj(intChar -> (char)intChar)
                .filter(base -> !matcher.reset(String.valueOf(base)).matches())
                .collect(Collectors.toSet());

        throw new InvalidBaseCharacterException(bases, invalidBaseChars);
    }

    private SAMRecord createSamRecord(boolean paired, String baseName, String read, String qualities) {

        final byte[] normalizedQualities = qualities.getBytes(StandardCharsets.UTF_8);
        qualityNormalizer.normalize(normalizedQualities);

        final SAMRecord rec = new SAMRecord(writer.getFileHeader());
        rec.setReadUnmappedFlag(true);
        rec.setAttribute(ReservedTagConstants.READ_GROUP_ID, DEFAULT_READ_GROUP_NAME);
        rec.setReadName(baseName);
        rec.setReadString(read);
        rec.setBaseQualities(normalizedQualities);

        if (paired) {
            rec.setReadPairedFlag(true);
            rec.setMateUnmappedFlag(true);
        }
        return rec;
    }
}
