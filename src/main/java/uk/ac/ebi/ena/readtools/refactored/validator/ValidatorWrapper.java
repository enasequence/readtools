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
package uk.ac.ebi.ena.readtools.refactored.validator;

import java.io.File;
import java.util.List;

import uk.ac.ebi.ena.readtools.refactored.FileFormat;
import uk.ac.ebi.ena.readtools.refactored.provider.FastqReadsProvider;
import uk.ac.ebi.ena.readtools.refactored.provider.ReadsProvider;
import uk.ac.ebi.ena.readtools.refactored.provider.SamReadsProvider;

public class ValidatorWrapper {
    protected final List<File> files;
    protected final FileFormat format;

    public ValidatorWrapper(List<File> files, FileFormat format) {
        this.files = files;
        this.format = format;
    }

    public void run() throws ReadsValidationException {
        switch (format) {
            case FASTQ:
                for (File file : files) {
                    validateFastq(file);
                }
            default:
                throw new ReadsValidationException("not implemented", 0);
        }
    }

    public static void validateFastq(File file) throws ReadsValidationException {
        try (ReadsProvider producer = new FastqReadsProvider(file)) {
            new InsdcReadsValidator().validate(producer);
            new FastqReadsValidator().validate(producer);
        } catch (ReadsValidationException rve) {
            throw rve;
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }

    public static void validateSam(File file) throws ReadsValidationException {
        try (ReadsProvider producer = new SamReadsProvider(file)) {
            new InsdcReadsValidator().validate(producer);
//            new SamReadsValidator().validate(producer);
        } catch (ReadsValidationException rve) {
            throw rve;
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }
}
