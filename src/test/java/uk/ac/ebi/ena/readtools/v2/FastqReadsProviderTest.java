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
package uk.ac.ebi.ena.readtools.v2;

import static org.junit.Assert.assertEquals;
import static uk.ac.ebi.ena.readtools.v2.TestFileUtil.createOutputFolder;
import static uk.ac.ebi.ena.readtools.v2.TestFileUtil.saveRandomized;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;

import org.junit.Test;

import htsjdk.samtools.util.FastqQualityFormat;

import uk.ac.ebi.ena.readtools.v2.provider.FastqReadsProvider;
import uk.ac.ebi.ena.readtools.v2.validator.ReadsValidationException;

public class FastqReadsProviderTest {
    @Test
    public void detectFormatNoNormalisation() throws IOException, ReadsValidationException {
        File output_dir = createOutputFolder();
        Path f1 = saveRandomized(
                "@R1\n"
                        + "ATCCATCC\n"
                        + "+\n"
                        + "!\"#${|}~\n",
                output_dir.toPath(), true, "fastq", "gz");

        FastqReadsProvider mrp = new FastqReadsProvider(f1.toFile(), false);
        assertEquals(null, mrp.getQualityFormat());
        assertEquals("!\"#${|}~", mrp.iterator().next().getQualityScores());
    }

    @Test
    public void detectFormatSolexa() throws IOException, ReadsValidationException {
        File output_dir = createOutputFolder();
        Path file = saveRandomized(
                "@R1\n"
                        + "ATCCATCC\n"
                        + "+\n"
                        + ";<=>{|}~\n",
                output_dir.toPath(), true, "fastq", "gz");

        FastqReadsProvider mrp = new FastqReadsProvider(file.toFile(), true);
        assertEquals(FastqQualityFormat.Solexa, mrp.getQualityFormat());
        assertEquals("!!!!\\]^_", mrp.iterator().next().getQualityScores());
    }

    @Test
    public void detectFormatIllumina() throws IOException, ReadsValidationException {
        File output_dir = createOutputFolder();
        Path file = saveRandomized(
                "@R1\n"
                        + "ATCCATCC\n"
                        + "+\n"
                        + "@ABC{|}~\n",
                output_dir.toPath(), true, "fastq", "gz");

        FastqReadsProvider mrp = new FastqReadsProvider(file.toFile(), true);
        assertEquals(FastqQualityFormat.Illumina, mrp.getQualityFormat());
        assertEquals("!\"#$\\]^_", mrp.iterator().next().getQualityScores());
    }

    @Test
    public void detectFormatStandard() throws IOException, ReadsValidationException {
        File output_dir = createOutputFolder();
        Path file = saveRandomized(
                "@R1\n"
                        + "ATCCATCC\n"
                        + "+\n"
                        + "!\"#${|}~\n",
                output_dir.toPath(), true, "fastq", "gz");

        FastqReadsProvider mrp = new FastqReadsProvider(file.toFile(), true);
        assertEquals(FastqQualityFormat.Standard, mrp.getQualityFormat());
        assertEquals("!\"#${|}~", mrp.iterator().next().getQualityScores());
    }
}
