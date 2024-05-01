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
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.stream.Collectors;
import uk.ac.ebi.ena.readtools.v2.FileFormat;
import uk.ac.ebi.ena.readtools.v2.validator.ReadsValidationException;
import uk.ac.ebi.ena.readtools.v2.validator.ValidatorWrapper;
import uk.ac.ebi.ena.readtools.webin.cli.rawreads.RawReadsFile;
import uk.ac.ebi.ena.readtools.webin.cli.rawreads.RawReadsFile.Filetype;
import uk.ac.ebi.ena.webin.cli.validator.api.ValidationResponse.status;
import uk.ac.ebi.ena.webin.cli.validator.api.Validator;
import uk.ac.ebi.ena.webin.cli.validator.file.SubmissionFile;
import uk.ac.ebi.ena.webin.cli.validator.file.SubmissionFiles;
import uk.ac.ebi.ena.webin.cli.validator.manifest.ReadsManifest;
import uk.ac.ebi.ena.webin.cli.validator.manifest.ReadsManifest.FileType;
import uk.ac.ebi.ena.webin.cli.validator.message.ValidationMessage;
import uk.ac.ebi.ena.webin.cli.validator.message.ValidationResult;
import uk.ac.ebi.ena.webin.cli.validator.response.ReadsValidationResponse;

public class ReadsValidator implements Validator<ReadsManifest, ReadsValidationResponse> {
  private static final long QUICK_READ_LIMIT = 100_000;
  private static final long EXTENDED_READ_LIMIT = 100_000_000;

  private List<ValidatorWrapper.FileQualityStats> fileQualityStats = new ArrayList<>();

  public List<ValidatorWrapper.FileQualityStats> getFileQualityStats() {
    return fileQualityStats;
  }

  @Override
  public ReadsValidationResponse validate(ReadsManifest manifest) {
    if (manifest == null) {
      throw new RuntimeException("Manifest is missing.");
    }
    if (manifest.getReportFile() == null) {
      throw new RuntimeException("Report file is missing.");
    }
    if (manifest.getProcessDir() == null) {
      throw new RuntimeException("Process directory is missing.");
    }

    return validate(
        manifest.getReportFile(), manifest.getQualityScore(), manifest.files(), manifest.isQuick());
  }

  public ReadsValidationResponse validate(
      File reportFile,
      ReadsManifest.QualityScore qualityScore,
      SubmissionFiles<FileType> submissionFiles,
      boolean isQuick) {
    ValidationResult result = new ValidationResult(reportFile);
    List<RawReadsFile> files = submissionFilesToRawReadsFiles(qualityScore, submissionFiles);

    return validate(result, files, isQuick);
  }

  public ReadsValidationResponse validate(
      ValidationResult result, List<RawReadsFile> files, boolean isQuick) {
    fileQualityStats.clear();

    AtomicBoolean paired = new AtomicBoolean();

    if (files != null && files.size() > 0) {
      Filetype fileType = files.get(0).getFiletype();

      if (Filetype.fastq.equals(fileType)
          || Filetype.bam.equals(fileType)
          || Filetype.cram.equals(fileType)) {
        runValidatorWrapper(
            result,
            files.stream().map(e -> new File(e.getFilename())).collect(Collectors.toList()),
            FileFormat.valueOf(fileType.toString().toUpperCase()),
            isQuick);
      } else {
        throw new RuntimeException("Unsupported file type: " + fileType.name());
      }
    }

    ReadsValidationResponse resp = new ReadsValidationResponse();
    resp.setStatus(result.isValid() ? status.VALIDATION_SUCCESS : status.VALIDATION_ERROR);
    resp.setPaired(paired.get());
    return resp;
  }

  private void runValidatorWrapper(
      ValidationResult result, List<File> files, FileFormat fileFormat, boolean isQuick) {
    ValidatorWrapper validatorWrapper =
        new ValidatorWrapper(files, fileFormat, isQuick ? QUICK_READ_LIMIT : EXTENDED_READ_LIMIT);
    try {
      validatorWrapper.run();
      fileQualityStats = validatorWrapper.getFileQualityStats();
    } catch (ReadsValidationException e) {
      result.add(ValidationMessage.error(e.getMessage()));
      e.printStackTrace();
    }
  }

  private List<RawReadsFile> submissionFilesToRawReadsFiles(
      ReadsManifest.QualityScore qualityScore, SubmissionFiles<FileType> submissionFiles) {
    List<RawReadsFile> files =
        submissionFiles.get().stream()
            .map(ReadsValidator::submissionFilesToRawReadsFile)
            .collect(Collectors.toList());

    RawReadsFile.AsciiOffset asciiOffset = null;
    RawReadsFile.QualityScoringSystem qualityScoringSystem = null;

    if (qualityScore != null) {
      switch (qualityScore) {
        case PHRED_33:
          asciiOffset = RawReadsFile.AsciiOffset.FROM33;
          qualityScoringSystem = RawReadsFile.QualityScoringSystem.phred;
          break;
        case PHRED_64:
          asciiOffset = RawReadsFile.AsciiOffset.FROM64;
          qualityScoringSystem = RawReadsFile.QualityScoringSystem.phred;
          break;
        case LOGODDS:
          //                    asciiOffset = null;
          qualityScoringSystem = RawReadsFile.QualityScoringSystem.log_odds;
          break;
      }
    }

    // Set FASTQ quality scoring system and ascii offset.

    for (RawReadsFile f : files) {
      if (f.getFiletype().equals(Filetype.fastq)) {
        if (qualityScoringSystem != null) {
          f.setQualityScoringSystem(qualityScoringSystem);
        }
        if (asciiOffset != null) {
          f.setAsciiOffset(asciiOffset);
        }
      }
    }

    return files;
  }

  public static RawReadsFile submissionFilesToRawReadsFile(SubmissionFile<FileType> file) {
    Path inputDir = file.getFile().toPath().getParent();

    RawReadsFile f = new RawReadsFile();
    f.setInputDir(inputDir);
    f.setReportFile(file.getReportFile().toPath());
    f.setFiletype(Filetype.valueOf(file.getFileType().name().toLowerCase()));

    String fileName = file.getFile().getPath();

    Path path = Paths.get(fileName);
    if (!path.isAbsolute()) {
      f.setFilename(inputDir.resolve(path).toString());
    } else {
      f.setFilename(fileName);
    }

    return f;
  }
}
