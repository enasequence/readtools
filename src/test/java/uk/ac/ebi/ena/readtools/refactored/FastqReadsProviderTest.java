package uk.ac.ebi.ena.readtools.refactored;

import htsjdk.samtools.util.FastqQualityFormat;
import org.junit.Test;
import uk.ac.ebi.ena.readtools.refactored.provider.FastqReadsProvider;
import uk.ac.ebi.ena.readtools.refactored.validator.ReadsValidationException;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;

import static org.junit.Assert.assertEquals;
import static uk.ac.ebi.ena.readtools.refactored.TestFileUtil.createOutputFolder;
import static uk.ac.ebi.ena.readtools.refactored.TestFileUtil.saveRandomized;

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

        FastqReadsProvider mrp = new FastqReadsProvider(f1.toFile());
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
