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
import java.util.ArrayList;
import java.util.List;

import uk.ac.ebi.ena.readtools.v2.FileFormat;
import uk.ac.ebi.ena.readtools.v2.provider.ReadsProviderFactory;

public class ValidatorWrapper {
    private final List<File> files;
    private final FileFormat format;
    private final long readCountLimit;

    private List<FileQualityStats> fileQualityStats = new ArrayList<>();

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

    public void run() throws ReadsValidationException {
        fileQualityStats.clear();

        switch (format) {
            case FASTQ:
                for (File file : files) {
                    validateFastq(file);
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
}
