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
package uk.ac.ebi.ena.readtools.validator;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;

import org.junit.Assert;
import org.junit.Test;

import uk.ac.ebi.ena.webin.cli.validator.api.ValidationResponse;
import uk.ac.ebi.ena.webin.cli.validator.file.SubmissionFile;
import uk.ac.ebi.ena.webin.cli.validator.file.SubmissionFiles;
import uk.ac.ebi.ena.webin.cli.validator.manifest.ReadsManifest;
import uk.ac.ebi.ena.webin.cli.validator.response.ReadsValidationResponse;

public class ReadsValidatorTest {

    @Test
    public void test() throws IOException {
        Path workDir = Files.createTempDirectory("rvt");
        Path reportFile = Files.createTempFile(workDir, "", ".report");

        File filePath1 = new File(ReadsValidatorTest.class.getClassLoader().getResource(
            "10x/3fastq/I1.fastq").getFile());
        File filePath2 = new File(ReadsValidatorTest.class.getClassLoader().getResource(
            "10x/3fastq/R1.fastq").getFile());
        File filePath3 = new File(ReadsValidatorTest.class.getClassLoader().getResource(
            "10x/3fastq/R2.fastq").getFile());

        SubmissionFiles<ReadsManifest.FileType> submissionFiles = new SubmissionFiles<>();
        submissionFiles.add(new SubmissionFile<>(ReadsManifest.FileType.FASTQ, filePath1, reportFile.toFile()));
        submissionFiles.add(new SubmissionFile<>(ReadsManifest.FileType.FASTQ, filePath2, reportFile.toFile()));
        submissionFiles.add(new SubmissionFile<>(ReadsManifest.FileType.FASTQ, filePath3, reportFile.toFile()));

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
