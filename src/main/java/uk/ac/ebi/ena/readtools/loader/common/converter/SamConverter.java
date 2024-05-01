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

import htsjdk.samtools.*;
import java.io.EOFException;
import java.io.IOException;
import java.io.InputStream;
import uk.ac.ebi.ena.readtools.cram.ref.ENAReferenceSource;
import uk.ac.ebi.ena.readtools.loader.common.writer.ReadWriter;
import uk.ac.ebi.ena.readtools.loader.fastq.Read;
import uk.ac.ebi.ena.readtools.loader.fastq.SamRecordWrapper;

public class SamConverter extends AbstractReadConverter<Read> {

  private final boolean isCram;
  private SAMRecordIterator samRecordIterator;

  public SamConverter(boolean isCram, InputStream istream, ReadWriter<Read, ?> writer) {
    super(istream, writer);
    this.isCram = isCram;
  }

  public SamConverter(
      boolean isCram, InputStream istream, ReadWriter<Read, ?> writer, Long readLimit) {
    super(istream, writer, readLimit);
    this.isCram = isCram;
  }

  @Override
  protected void begin() {
    try {
      SamInputResource samInputResource = SamInputResource.of(istream);

      SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);
      SamReaderFactory factory = SamReaderFactory.make();
      factory.enable(SamReaderFactory.Option.DONT_MEMORY_MAP_INDEX);
      factory.validationStringency(ValidationStringency.SILENT);
      factory.samRecordFactory(DefaultSAMRecordFactory.getInstance());

      if (isCram) {
        ENAReferenceSource reference_source = new ENAReferenceSource();
        factory.referenceSource(reference_source);
      }
      SamReader samReader = factory.open(samInputResource);
      samRecordIterator = samReader.iterator();
    } catch (Exception ex) {
      throw new ConverterException(ex);
    }
  }

  @Override
  public Read getNextSpotFromInputStream(InputStream inputStream) throws IOException {
    if (samRecordIterator.hasNext()) {
      SAMRecord samRecord = samRecordIterator.next();
      return new SamRecordWrapper(samRecord);
    } else {
      throw new EOFException();
    }
  }
}
