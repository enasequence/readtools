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

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Paths;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import htsjdk.samtools.ReservedTagConstants;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;

import uk.ac.ebi.ena.readtools.common.reads.QualityNormalizer;
import uk.ac.ebi.ena.readtools.loader.common.InvalidBaseCharacterException;
import uk.ac.ebi.ena.readtools.loader.common.writer.ReadWriter;
import uk.ac.ebi.ena.readtools.loader.common.writer.ReadWriterException;
import uk.ac.ebi.ena.readtools.loader.common.writer.Spot;
import uk.ac.ebi.ena.readtools.loader.fastq.PairedRead;
import uk.ac.ebi.ena.readtools.loader.fastq.Read;
import uk.ac.ebi.ena.readtools.utils.Utils;

/**
 * Accepts Fastq spot data and writes them out to a BAM file.
 */
public class Fastq2BamWriter implements ReadWriter<PairedRead, Spot> {
    private static final String DEFAULT_READ_GROUP_NAME = "A";
    private static final String VALID_DNA_CHARSET = ".acmgrsvtwyhkdbnACMGRSVTWYHKDBN";

    private final QualityNormalizer qualityNormalizer;
    private final String sampleName;
    private final boolean convertUracil;
    private final boolean paired;

    private final Pattern validDnaCharsetPattern;
    private final SAMFileWriter writer;

    public Fastq2BamWriter(QualityNormalizer qualityNormalizer, String sampleName, String outputFilePath,
                           String tempDir, boolean convertUracil, boolean paired) {
        this.qualityNormalizer = qualityNormalizer;
        this.sampleName = sampleName;
        this.convertUracil = convertUracil;
        this.paired = paired;

        if (sampleName == null || sampleName.trim().isEmpty()) {
            throw new IllegalArgumentException("Sample name is either null or empty.");
        }

        validDnaCharsetPattern = convertUracil ? Pattern.compile("^[" + VALID_DNA_CHARSET + "uU]+$")
                : Pattern.compile("^[" + VALID_DNA_CHARSET + "]+$");

        writer = new SAMFileWriterFactory().setTempDirectory(new File(tempDir)).makeSAMOrBAMWriter(
                createHeader(), false, Paths.get(outputFilePath));
    }

    @Override
    public void write(PairedRead spot) throws ReadWriterException {
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
                Read unpaired = spot.getUnpaired();
                SAMRecord rec = createSamRecord(false, spot.name, unpaired.bases, unpaired.quals);
                rec.setReadPairedFlag(false);
                writer.addAlignment(rec);
            }
        } catch (Exception ex) {
            throw new ReadWriterException(ex);
        }
    }

    @Override
    public void cascadeErrors() throws ReadWriterException {
    }

    @Override
    public void setWriter(ReadWriter<Spot, ? extends Spot> readWriter) {
    }

    public void unwind() {
        writer.close();
    }

    private SAMFileHeader createHeader() {
        final SAMReadGroupRecord rgroup = new SAMReadGroupRecord(DEFAULT_READ_GROUP_NAME);
        rgroup.setSample(this.sampleName);

        final SAMFileHeader header = new SAMFileHeader();
        header.addReadGroup(rgroup);
        if (paired) {
            header.setSortOrder(SAMFileHeader.SortOrder.queryname);
        } else {
            header.setSortOrder(SAMFileHeader.SortOrder.unsorted);
        }

        return header;
    }

    private void validate(PairedRead pairedRead) {
        if (pairedRead.forward != null) {
            Matcher matcher = validDnaCharsetPattern.matcher(pairedRead.forward.bases);
            if (!matcher.matches()) {
                handleInvalidDnaCharset(pairedRead.forward.bases, matcher);
            }

            if (pairedRead.forward.bases.length() != pairedRead.forward.quals.length())
                throw new IllegalArgumentException(String.format(
                        "FATAL: Spot bases and qualities length do not match. Malformed spot\n%s\n", pairedRead));
        }

        if (pairedRead.reverse != null) {
            Matcher matcher = validDnaCharsetPattern.matcher(pairedRead.reverse.bases);
            if (!matcher.matches())
                handleInvalidDnaCharset(pairedRead.reverse.bases, matcher);

            if (pairedRead.reverse.bases.length() != pairedRead.reverse.quals.length())
                throw new IllegalArgumentException(String.format(
                        "FATAL: Spot bases and qualities length do not match. Malformed spot\n%s\n", pairedRead));
        }
    }

    private void handleInvalidDnaCharset(String bases, Matcher matcher) {
        Set<Character> invalidBaseChars = bases.chars()
                .mapToObj(intChar -> (char) intChar)
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
        rec.setReadString(modifyBases(read));
        rec.setBaseQualities(normalizedQualities);

        if (paired) {
            rec.setReadPairedFlag(true);
            rec.setMateUnmappedFlag(true);
        }
        return rec;
    }

    private String modifyBases(String bases) {
        return convertUracil ? Utils.replaceUracilBases(bases) : bases;
    }
}
