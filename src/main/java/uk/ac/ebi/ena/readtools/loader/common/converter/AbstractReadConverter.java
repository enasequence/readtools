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

import java.io.BufferedInputStream;
import java.io.EOFException;
import java.io.IOException;
import java.io.InputStream;
import uk.ac.ebi.ena.readtools.loader.common.writer.ReadWriter;
import uk.ac.ebi.ena.readtools.loader.common.writer.ReadWriterException;
import uk.ac.ebi.ena.readtools.loader.common.writer.ReadWriterMemoryLimitException;
import uk.ac.ebi.ena.readtools.loader.common.writer.Spot;

public abstract class AbstractReadConverter<T extends Spot> implements Converter {
  protected final InputStream istream;
  protected final Long readLimit;
  protected final ReadWriter<T, ?> readWriter;

  long readCount = 0, baseCount = 0;

  protected boolean isEofReached = false;

  protected AbstractReadConverter(InputStream istream, ReadWriter<T, ?> writer) {
    this(istream, writer, null);
  }

  /**
   * @param istream
   * @param readLimit Only read limited amount of reads.
   */
  protected AbstractReadConverter(InputStream istream, ReadWriter<T, ?> writer, Long readLimit) {
    this.istream = new BufferedInputStream(istream, 1024 * 1024);
    this.readWriter = writer;
    this.readLimit = readLimit;
  }

  /**
   * Get the total number of reads that were read.
   *
   * @return
   */
  public long getReadCount() {
    return readCount;
  }

  /**
   * Get the total number of bases that were read.
   *
   * @return
   */
  public long getBaseCount() {
    return baseCount;
  }

  public final void run() {
    try {
      begin();

      do {
        readWriter.write(getNextSpot());
      } while (!isDone());
    } catch (ConverterEOFException ignored) {
      isEofReached = true;
    } catch (Exception e) {
      if (e instanceof ReadWriterException
          || e instanceof ReadWriterMemoryLimitException
          || e instanceof ConverterException) {
        throw e;
      } else {
        throw new RuntimeException(e);
      }
    } finally {
      end();
    }
  }

  public void runOnce() {
    try {
      if (!isDone()) {
        readWriter.write(getNextSpot());
      }
    } catch (ConverterEOFException ignored) {
      isEofReached = true;
    } catch (Exception e) {
      if (e instanceof ReadWriterException
          || e instanceof ReadWriterMemoryLimitException
          || e instanceof ConverterException) {
        throw e;
      } else {
        throw new RuntimeException(e);
      }
    }
  }

  public boolean isDone() {
    return !isWithinReadLimit() || isEofReached;
  }

  protected void begin() {}

  protected void end() {}

  private boolean isWithinReadLimit() {
    if (readLimit == null) {
      return true;
    } else {
      return readCount < readLimit;
    }
  }

  // Re-implement if you need special type of feeding
  private T getNextSpot() {
    T spot;

    try {
      spot = getNextSpotFromInputStream(istream);
      ++readCount;
      baseCount += spot.getBaseCount();

      return spot;
    } catch (EOFException e) {
      throw new ConverterEOFException(readCount);
    } catch (ConverterException e) {
      throw e;
    } catch (Throwable cause) {
      throw new ConverterException(cause);
    }
  }

  public abstract T getNextSpotFromInputStream(InputStream inputStream) throws IOException;
}
