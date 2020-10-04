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
import uk.ac.ebi.ena.readtools.loader.common.consumer.DataConsumable;
import uk.ac.ebi.ena.readtools.loader.common.consumer.DataConsumer;
import uk.ac.ebi.ena.readtools.loader.common.consumer.DataConsumerException;
import uk.ac.ebi.ena.readtools.loader.fastq.IlluminaSpot;

import java.nio.charset.StandardCharsets;
import java.nio.file.Paths;
import java.util.Set;
import java.util.function.Consumer;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

/**
 * Accepts Illumina spot data and writes them out to a BAM file.
 */
public class Fastq2BamConsumer implements DataConsumer<IlluminaSpot, DataConsumable> {

    private static final String DEFAULT_READ_GROUP_NAME = "A";

    private static final String VALID_DNA_CHARSET = "^[.acmgrsvtwyhkdbnNACMGRSVTWYHKDBN]+$";
    private static final Pattern VALID_DNA_CHARSET_PATTERN = Pattern.compile(VALID_DNA_CHARSET);

    private final QualityNormalizer qualityNormalizer;

    private final String sampleName;

    private final SAMFileWriter writer;

    private volatile boolean isOk = true;

    public Fastq2BamConsumer(QualityNormalizer qualityNormalizer, String sampleName, String outputFilePath) {
        this.qualityNormalizer = qualityNormalizer;
        this.sampleName = sampleName;

        if (sampleName == null || sampleName.trim().isEmpty()) {
            throw new IllegalArgumentException("Sample name is either null or empty.");
        }

        writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(
                createHeader(), false, Paths.get(outputFilePath));
    }

    @Override
    public void consume(IlluminaSpot iSpot) throws DataConsumerException {
        try {
            Matcher matcher = VALID_DNA_CHARSET_PATTERN.matcher(iSpot.bases);
            if( !matcher.matches() )
                handleInvalidDnaCharset(iSpot.bases, matcher);

            if( iSpot.bases.length() != iSpot.quals.length() )
                throw new IllegalArgumentException( String.format( "FATAL: Spot bases and qualities length do not match. Malformed spot\n%s\n", iSpot ) );

            if (isPaired(iSpot)) {
                SAMRecord rec1 = createSamRecord(
                        true, iSpot.name, getForwardBases(iSpot), getForwardQualities(iSpot));
                rec1.setFirstOfPairFlag(true);
                rec1.setSecondOfPairFlag(false);
                writer.addAlignment(rec1);

                SAMRecord rec2 = createSamRecord(
                        true, iSpot.name, getReverseBases(iSpot), getReverseQualities(iSpot));
                rec2.setFirstOfPairFlag(false);
                rec2.setSecondOfPairFlag(true);
                writer.addAlignment(rec2);

            } else {
                SAMRecord rec = createSamRecord(false, iSpot.name, iSpot.bases, iSpot.quals);
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
    public void setConsumer(DataConsumer<DataConsumable, ? extends DataConsumable> dataConsumer) { }

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

    private void handleInvalidDnaCharset(String bases, Matcher matcher) {
        Set<Character> invalidBaseChars = bases.chars()
                .mapToObj(intChar -> (char)intChar)
                .filter(base -> !matcher.reset(String.valueOf(base)).matches())
                .collect(Collectors.toSet());

        throw new InvalidBaseCharacterException(bases, invalidBaseChars);
    }

    private boolean isPaired(IlluminaSpot iSpot) {
        return iSpot.read_name.length == 2
                && iSpot.read_name[IlluminaSpot.FORWARD] != null
                && iSpot.read_name[IlluminaSpot.REVERSE] != null;
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

    private String getForwardBases(IlluminaSpot iSpot) {
        int start = iSpot.read_start[IlluminaSpot.FORWARD];
        int end = iSpot.read_length[IlluminaSpot.FORWARD];

        return iSpot.bases.substring(start, end);
    }

    private String getForwardQualities(IlluminaSpot iSpot) {
        int start = iSpot.read_start[IlluminaSpot.FORWARD];
        int end = iSpot.read_length[IlluminaSpot.FORWARD];

        return iSpot.quals.substring(start, end);
    }

    private String getReverseBases(IlluminaSpot iSpot) {
        int start = iSpot.read_start[IlluminaSpot.REVERSE];

        return iSpot.bases.substring(start);
    }

    private String getReverseQualities(IlluminaSpot iSpot) {
        int start = iSpot.read_start[IlluminaSpot.REVERSE];

        return iSpot.quals.substring(start);
    }
}
