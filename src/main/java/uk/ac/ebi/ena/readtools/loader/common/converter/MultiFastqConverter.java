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

import java.io.EOFException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.List;
import uk.ac.ebi.ena.readtools.common.reads.QualityNormalizer;
import uk.ac.ebi.ena.readtools.loader.common.writer.ReadWriter;
import uk.ac.ebi.ena.readtools.loader.common.writer.ReadWriterException;
import uk.ac.ebi.ena.readtools.loader.common.writer.ReadWriterMemoryLimitException;
import uk.ac.ebi.ena.readtools.loader.common.writer.Spot;
import uk.ac.ebi.ena.readtools.loader.fastq.Read;

/**
 * Similar to {@link AutoNormalizeQualityReadConverter}, but here, base quality normalizer is
 * provided explicitly. If no normalizer is provided then normalization is not performed.
 */
public class MultiFastqConverter<T extends Spot> implements Converter {
  List<ReadReader> readers = new ArrayList<>();
  List<InputStream> istreams = new ArrayList<>();
  final ReadWriter<Read, T> readWriter;
  final Long readLimit;

  long readCount = 0, baseCount = 0;

  List<ReadReader> readersCompleted = new ArrayList<>();
  List<InputStream> istreamsCompleted = new ArrayList<>();

  int readerIndex;

  public MultiFastqConverter(List<InputStream> istreams, ReadWriter<Read, T> readWriter) {
    this(istreams, readWriter, null);
  }

  public MultiFastqConverter(
      List<InputStream> istreams, ReadWriter<Read, T> readWriter, Long readLimit) {
    this.readWriter = readWriter;

    for (int readerIndex = 0; readerIndex < istreams.size(); readerIndex++) {
      this.istreams.add(istreams.get(readerIndex));
      readers.add(new ReadReader(String.valueOf(readerIndex + 1)));
    }

    this.readLimit = readLimit;
  }

  public MultiFastqConverter(
      List<InputStream> istreams,
      List<QualityNormalizer> normalizers,
      ReadWriter<Read, T> readWriter,
      Long readLimit) {

    this.readWriter = readWriter;

    for (int readerIndex = 0; readerIndex < istreams.size(); readerIndex++) {
      this.istreams.add(istreams.get(readerIndex));
      readers.add(
          new ReadReader(
              normalizers.get(normalizers.size() == istreams.size() ? readerIndex : 0),
              String.valueOf(readerIndex + 1)));
    }

    this.readLimit = readLimit;
  }

  private boolean isWithinReadLimit() {
    if (readLimit == null) {
      return true;
    } else {
      return readCount < readLimit;
    }
  }

  public long getReadCount() {
    return readCount;
  }

  public long getBaseCount() {
    return baseCount;
  }

  public void run() {
    try {
      do {
        readWriter.write(convert());
      } while (!isDone());
    } catch (ReadWriterException e) {
      throw new ConverterPanicException(e);
    } catch (ConverterEOFException ignored) {
    } catch (Exception e) {
      if (e instanceof ReadWriterMemoryLimitException || e instanceof ConverterException) {
        throw e;
      } else {
        throw new RuntimeException(e);
      }
    }
  }

  public void runOnce() {
    try {
      if (!isDone()) {
        readWriter.write(convert());
      }
    } catch (ReadWriterException e) {
      throw new ConverterPanicException(e);
    } catch (ConverterEOFException ignored) {
    } catch (Exception e) {
      if (e instanceof ReadWriterMemoryLimitException || e instanceof ConverterException) {
        throw e;
      } else {
        throw new RuntimeException(e);
      }
    }
  }

  public boolean isDone() {
    return (readers.isEmpty() || !isWithinReadLimit());
  }

  private Read convert() {
    try {
      if (readerIndex >= readers.size()) {
        readerIndex = 0;
      }

      Read spot = readers.get(readerIndex).read(istreams.get(readerIndex));

      readerIndex++;

      readCount++;
      baseCount += spot.getBaseCount();

      return spot;
    } catch (EOFException e) {
      readersCompleted.add(readers.get(readerIndex));
      readers.remove(readerIndex);

      istreamsCompleted.add(istreams.get(readerIndex));
      istreams.remove(readerIndex);

      throw new ConverterEOFException(readCount);
    } catch (ConverterException e) {
      throw e;
    } catch (Exception cause) {
      throw new ConverterException(cause);
    }
  }
}
