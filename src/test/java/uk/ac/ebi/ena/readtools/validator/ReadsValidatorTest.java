package uk.ac.ebi.ena.readtools.validator;

import org.junit.Assert;
import org.junit.Test;
import uk.ac.ebi.ena.webin.cli.validator.api.ValidationResponse;
import uk.ac.ebi.ena.webin.cli.validator.file.SubmissionFile;
import uk.ac.ebi.ena.webin.cli.validator.file.SubmissionFiles;
import uk.ac.ebi.ena.webin.cli.validator.manifest.ReadsManifest;
import uk.ac.ebi.ena.webin.cli.validator.response.ReadsValidationResponse;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;

public class ReadsValidatorTest {

    @Test
    public void test() throws IOException {
        Path workDir = Files.createTempDirectory("rvt");
        Path reportFile = Files.createTempFile(workDir, "", ".report");

        File filePath1 = new File(ReadsValidatorTest.class.getClassLoader().getResource(
            "uk/ac/ebi/ena/webin/cli/rawreads/EP0_GTTCCTT_S1.txt.gz").getFile());
        File filePath2 = new File(ReadsValidatorTest.class.getClassLoader().getResource(
            "uk/ac/ebi/ena/webin/cli/rawreads/EP0_GTTCCTT_S2.txt.gz").getFile());

        SubmissionFiles<ReadsManifest.FileType> submissionFiles = new SubmissionFiles<>();
        submissionFiles.add(new SubmissionFile<>(ReadsManifest.FileType.FASTQ, filePath1, reportFile.toFile()));
        submissionFiles.add(new SubmissionFile<>(ReadsManifest.FileType.FASTQ, filePath2, reportFile.toFile()));

        ReadsManifest readsManifest = new ReadsManifest();
        readsManifest.setReportFile(reportFile.toFile());
        readsManifest.setProcessDir(workDir.toFile());
        readsManifest.setFiles(submissionFiles);
        readsManifest.setQuick(true);

        ReadsValidator readsValidator = new ReadsValidator();

        ReadsValidationResponse validationResponse = readsValidator.validate(readsManifest);

        Assert.assertEquals(ValidationResponse.status.VALIDATION_SUCCESS, validationResponse.getStatus());
    }
}
