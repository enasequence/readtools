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
package uk.ac.ebi.ena.readtools.v2.validator;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;

import uk.ac.ebi.ena.readtools.v2.FileFormat;
import uk.ac.ebi.ena.readtools.v2.provider.ReadsProviderFactory;
import uk.ac.ebi.ena.readtools.webin.cli.rawreads.BloomWrapper;
import uk.ac.ebi.ena.webin.cli.validator.message.ValidationMessage;

import static uk.ac.ebi.ena.readtools.v2.validator.FastqReadsValidator.EXPECTED_SIZE;

public class ValidatorWrapper {
    private final List<File> files;
    private final FileFormat format;
    private final long readCountLimit;
    private static final int PAIRING_THRESHOLD = 20;

    private List<FileQualityStats> fileQualityStats = new ArrayList<>();
    private Set<String> labels = new HashSet<>();
    boolean paired = false;

    public static class FileQualityStats {
        private File file = null;
        private long readCount = 0;
        private long highQualityReadCount = 0;

        public FileQualityStats() {}

        public FileQualityStats(File file, long readCount, long highQualityReadCount) {
            this.file = file;
            this.readCount = readCount;
            this.highQualityReadCount = highQualityReadCount;
        }

        public File getFile() {
            return file;
        }

        public long getReadCount() {
            return readCount;
        }

        public long getHighQualityReadCount() {
            return highQualityReadCount;
        }
    }

    public ValidatorWrapper(List<File> files, FileFormat format, long readCountLimit) {
        this.files = files;
        this.format = format;
        this.readCountLimit = readCountLimit;
    }

    public List<FileQualityStats> getFileQualityStats() {
        return fileQualityStats;
    }

    public boolean isPaired() {
        return paired;
    }

    public void run() throws ReadsValidationException {
        fileQualityStats.clear();

        switch (format) {
            case FASTQ:
                if (files.size() == 1) {
                    validateFastq(files.get(0));
                } else {
                    validateMultipleFastq(files);
                }
                break;
            case BAM:
            case CRAM:
                for (File file : files) {
                    validateSam(file);
                }
                break;
            default:
                throw new ReadsValidationException("not implemented");
        }
    }

    public void validateFastq(File file) throws ReadsValidationException {
        try {
            ReadsProviderFactory factory = new ReadsProviderFactory(file, format);
            validateInsdc(file);
            new FastqReadsValidator(readCountLimit).validate(factory);
        } catch (ReadsValidationException rve) {
            throw rve;
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }

    private void validateMultipleFastq(List<File> files) throws ReadsValidationException {
        try {
            /** Should ideally have high number of duplicates as it will point to higher pairing percentage.
             * To keep memory consumption lower, because we can tolerate false positive here, use lower expected read size. */
            BloomWrapper mainFileOnlyPairingBloomWrapper = new BloomWrapper(EXPECTED_SIZE / 10);
            List<PairedFiles> pairedFiles = new ArrayList<>();
            int fileNumber = 0;
            for (File file : files) {
                ReadsProviderFactory factory = new ReadsProviderFactory(file, format);
                validateInsdc(file);

                BloomWrapper bloomWrapper;
                if (fileNumber == 0) {
                    bloomWrapper = mainFileOnlyPairingBloomWrapper;
                } else {
                    //Make a copy of the main file's pairing information so we do not have re-create it for every other file.
                    bloomWrapper = mainFileOnlyPairingBloomWrapper.getCopy();
                }
                FastqReadsValidator validator = new FastqReadsValidator(readCountLimit, bloomWrapper, labels);
                validator.validate(factory);
                fileNumber++;

                long readCount = Math.max(mainFileOnlyPairingBloomWrapper.getAddCount(), bloomWrapper.getAddCount());
                long pairedCount = bloomWrapper.getPossibleDuplicateCount();
                double pairingPercentage = 100 * ((double) pairedCount / (double) readCount);

                pairedFiles.add(new PairedFiles(
                        files.get(0).getAbsolutePath(), file.getAbsolutePath(), pairingPercentage));
            }

            paired = false;
            //Label set size and low pairing percentage validation.
            if (labels.size() <= files.size()) {
                paired = true;

                PairedFiles lowestPairingPercentagePair = pairedFiles.stream()
                        .sorted(Comparator.comparingDouble(pair -> pair.pairingPercentage))
                        .findFirst().get();

                //TODO: estimate bloom false positives impact
                if (lowestPairingPercentagePair.pairingPercentage < (double) PAIRING_THRESHOLD) {
                    //terminal error
//                    validationResult.add(ValidationMessage.error(
//                            String.format("Detected paired fastq submission with less than %d%% of paired reads between %s and %s",
//                                    PAIRING_THRESHOLD, lowestPairingPercentagePair.fileName1, lowestPairingPercentagePair.fileName2)));
                }
            } else {
//                validationResult.add(ValidationMessage.error(String.format(
//                        "When submitting paired reads using two Fastq files the reads must follow Illumina paired read naming conventions. "
//                                + "This was not the case for the submitted Fastq files: %s. Unable to determine pairing from set: %s",
//                        files,
//                        labels.stream().limit(10).collect(Collectors.joining(",", "", 10 < labels.size() ? "..." : "")))));
            }
        } catch (ReadsValidationException rve) {
            throw rve;
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }

    public void validateSam(File file) throws ReadsValidationException {
        try {
            validateInsdc(file);
        } catch (ReadsValidationException rve) {
            throw rve;
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }

    private void validateInsdc(File file) throws ReadsValidationException {
        ReadsProviderFactory factory = new ReadsProviderFactory(file, format);
        InsdcReadsValidator insdcReadsValidator = new InsdcReadsValidator(readCountLimit);
        insdcReadsValidator.validate(factory);

        fileQualityStats.add(new FileQualityStats(
                file, insdcReadsValidator.getReadCount(), insdcReadsValidator.getHighQualityReadCount()));
    }

    private static class PairedFiles {
        public String fileName1;
        public String fileName2;

        public double pairingPercentage;

        public PairedFiles(String fileName1, String fileName2, double pairingPercentage) {
            this.fileName1 = fileName1;
            this.fileName2 = fileName2;
            this.pairingPercentage = pairingPercentage;
        }
    }
}
