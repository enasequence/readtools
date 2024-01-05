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
import java.util.List;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.stream.Collectors;

import uk.ac.ebi.ena.readtools.refactored.FileFormat;
import uk.ac.ebi.ena.readtools.refactored.validator.ReadsValidationException;
import uk.ac.ebi.ena.readtools.refactored.validator.ValidatorWrapper;
import uk.ac.ebi.ena.readtools.webin.cli.rawreads.FastqScanner;
import uk.ac.ebi.ena.readtools.webin.cli.rawreads.RawReadsFile;
import uk.ac.ebi.ena.readtools.webin.cli.rawreads.RawReadsFile.Filetype;
import uk.ac.ebi.ena.webin.cli.validator.api.ValidationResponse.status;
import uk.ac.ebi.ena.webin.cli.validator.api.Validator;
import uk.ac.ebi.ena.webin.cli.validator.file.SubmissionFile;
import uk.ac.ebi.ena.webin.cli.validator.manifest.ReadsManifest;
import uk.ac.ebi.ena.webin.cli.validator.manifest.ReadsManifest.FileType;
import uk.ac.ebi.ena.webin.cli.validator.message.ValidationMessage;
import uk.ac.ebi.ena.webin.cli.validator.message.ValidationResult;
import uk.ac.ebi.ena.webin.cli.validator.response.ReadsValidationResponse;

public class ReadsValidator
        implements Validator<ReadsManifest, ReadsValidationResponse> {
    private static final long QUICK_READ_LIMIT = 100_000;
    private static final long EXTENDED_READ_LIMIT = 100_000_000;

    private ReadsManifest manifest;

    @Override
    public ReadsValidationResponse
    validate(ReadsManifest manifest) {
        this.manifest = manifest;
        if (manifest == null) {
            throw new RuntimeException("Manifest is missing.");
        }
        if (manifest.getReportFile() == null) {
            throw new RuntimeException("Report file is missing.");
        }
        if (manifest.getProcessDir() == null) {
            throw new RuntimeException("Process directory is missing.");
        }
        return validateSubmissionForContext();
    }

    private ReadsValidationResponse
    validateSubmissionForContext() {
        ValidationResult result = new ValidationResult(manifest.getReportFile());

        AtomicBoolean paired = new AtomicBoolean();

        List<RawReadsFile> files = createReadFiles();
        if (files != null && files.size() > 0) {
            Filetype fileType = files.get(0).getFiletype();

            if (files.size() == 2 && Filetype.fastq.equals(fileType)) {
                readFastqFile(result, files, paired);
            } else if (Filetype.fastq.equals(fileType)
                    || Filetype.bam.equals(fileType)
                    || Filetype.cram.equals(fileType)) {
                runValidatorWrapper(
                        result,
                        files.stream().map(e -> new File(e.getFilename())).collect(Collectors.toList()),
                        FileFormat.valueOf(fileType.toString().toUpperCase()));
            } else {
                throw new RuntimeException("Unsupported file type: " + fileType.name());
            }
        }

        ReadsValidationResponse resp = new ReadsValidationResponse();
        resp.setStatus(result.isValid() ? status.VALIDATION_SUCCESS : status.VALIDATION_ERROR);
        resp.setPaired(paired.get());
        return resp;
    }

    private void
    runValidatorWrapper(ValidationResult result, List<File> files, FileFormat fileFormat) {
        ValidatorWrapper validatorWrapper = new ValidatorWrapper(
                files, fileFormat, manifest.isQuick() ? QUICK_READ_LIMIT : EXTENDED_READ_LIMIT);
        try {
            validatorWrapper.run();
        } catch (ReadsValidationException e) {
            result.add(ValidationMessage.error(e.getMessage()));
            e.printStackTrace();
        }
    }

    private void
    readFastqFile(ValidationResult result, List<RawReadsFile> files, AtomicBoolean paired) {
        try {
            FastqScanner fs = new FastqScanner(manifest.isQuick() ? QUICK_READ_LIMIT : EXTENDED_READ_LIMIT) {
                @Override
                protected void logFlushMsg(String msg) {
                }

                @Override
                protected void logProcessedReadNumber(Long count) {
                    ;
                }
            };

            fs.checkFiles(result, files.toArray(new RawReadsFile[files.size()]));
            paired.set(fs.getPaired());
        } catch (Exception ex) {
            throw new RuntimeException(ex);
        } catch (Throwable throwable) {
            throw new RuntimeException(throwable);
        }
    }

    private List<RawReadsFile>
    createReadFiles() {
        RawReadsFile.AsciiOffset asciiOffset = null;
        RawReadsFile.QualityScoringSystem qualityScoringSystem = null;

        if (manifest.getQualityScore() != null) {
            switch (manifest.getQualityScore()) {
                case PHRED_33:
                    asciiOffset = RawReadsFile.AsciiOffset.FROM33;
                    qualityScoringSystem = RawReadsFile.QualityScoringSystem.phred;
                    break;
                case PHRED_64:
                    asciiOffset = RawReadsFile.AsciiOffset.FROM64;
                    qualityScoringSystem = RawReadsFile.QualityScoringSystem.phred;
                    break;
                case LOGODDS:
                    asciiOffset = null;
                    qualityScoringSystem = RawReadsFile.QualityScoringSystem.log_odds;
                    break;
            }
        }
        List<RawReadsFile> files = this.manifest.files().get()
                .stream()
                .map(file -> createRawReadsFile(file))
                .collect(Collectors.toList());

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


    public static RawReadsFile
    createRawReadsFile(SubmissionFile<FileType> file) {
        Path inputDir = file.getFile().toPath().getParent();

        RawReadsFile f = new RawReadsFile();
        f.setInputDir(inputDir);
        f.setReportFile(file.getReportFile().toPath());
        f.setFiletype(Filetype.valueOf(file.getFileType().name().toLowerCase()));

        String fileName = file.getFile().getPath();

        if (!Paths.get(fileName).isAbsolute()) {
            f.setFilename(inputDir.resolve(Paths.get(fileName)).toString());
        } else {
            f.setFilename(fileName);
        }

        return f;
    }
}
