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
package uk.ac.ebi.ena.readtools.loader.common.converter;

import htsjdk.samtools.util.FastqQualityFormat;
import java.io.IOException;
import java.io.InputStream;
import uk.ac.ebi.ena.readtools.common.reads.QualityNormalizer;
import uk.ac.ebi.ena.readtools.loader.common.writer.ReadWriter;
import uk.ac.ebi.ena.readtools.loader.fastq.Read;
import uk.ac.ebi.ena.readtools.utils.Utils;

/**
 * Similar to {@link ReadConverter}. except here, quality normalization strategy is chosen
 * automatically and quality normalization is always performed.
 */
public class AutoNormalizeQualityReadConverter extends AbstractReadConverter<Read> {

  private final String filePath;

  private final String defaultReadIndex;

  private volatile ReadReader readReader;

  public AutoNormalizeQualityReadConverter(
      InputStream istream, ReadWriter<Read, ?> writer, String defaultReadIndex, String filePath) {
    super(istream, writer);

    this.defaultReadIndex = defaultReadIndex;
    this.filePath = filePath;
  }

  /**
   * @param istream
   * @param readLimit Only read limited amount of reads.
   * @param defaultReadIndex
   * @param filePath
   */
  public AutoNormalizeQualityReadConverter(
      InputStream istream,
      ReadWriter<Read, ?> writer,
      Long readLimit,
      String defaultReadIndex,
      String filePath) {
    super(istream, writer, readLimit);

    this.defaultReadIndex = defaultReadIndex;
    this.filePath = filePath;
  }

  @Override
  protected void begin() {
    try {
      FastqQualityFormat qualityFormat = Utils.detectFastqQualityFormat(filePath, null);

      QualityNormalizer normalizer = Utils.getQualityNormalizer(qualityFormat);

      readReader = new ReadReader(normalizer, defaultReadIndex);
    } catch (Exception ex) {
      throw new ConverterException(ex);
    }
  }

  @Override
  public Read getNextSpotFromInputStream(InputStream inputStream) throws IOException {
    return readReader.read(inputStream);
  }
}
