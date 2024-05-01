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
package uk.ac.ebi.ena.readtools.v2.provider;

import java.io.File;
import uk.ac.ebi.ena.readtools.v2.FileFormat;
import uk.ac.ebi.ena.readtools.v2.read.IRead;
import uk.ac.ebi.ena.readtools.v2.validator.ReadsValidationException;

public class ReadsProviderFactory {
  private final File file;
  private final FileFormat format;
  private final boolean normaliseFastqQualityScores;

  public ReadsProviderFactory(File file, FileFormat format) {
    this(file, format, true);
  }

  public ReadsProviderFactory(File file, FileFormat format, boolean normaliseFastqQualityScores) {
    this.file = file;
    this.format = format;
    this.normaliseFastqQualityScores = normaliseFastqQualityScores;
  }

  public ReadsProvider<? extends IRead> makeReadsProvider() throws ReadsValidationException {
    switch (format) {
      case FASTQ:
        return new FastqReadsProvider(file, normaliseFastqQualityScores);
      case BAM:
      case CRAM:
        return new SamReadsProvider(file);
      default:
        throw new ReadsValidationException("not implemented");
    }
  }
}
