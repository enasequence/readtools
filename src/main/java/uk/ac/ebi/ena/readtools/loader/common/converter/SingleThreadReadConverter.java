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

import uk.ac.ebi.ena.readtools.common.reads.QualityNormalizer;
import uk.ac.ebi.ena.readtools.loader.common.writer.ReadWriter;
import uk.ac.ebi.ena.readtools.loader.common.writer.ReadWriterMemoryLimitException;
import uk.ac.ebi.ena.readtools.loader.common.writer.Spot;
import uk.ac.ebi.ena.readtools.loader.fastq.PairedFastqWriter;
import uk.ac.ebi.ena.readtools.loader.fastq.Read;

import java.io.EOFException;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.List;
import java.util.function.Supplier;

/**
 * Similar to {@link AutoNormalizeQualityReadConverter}, but here, base quality normalizer is provided explicitly. If no normalizer
 * is provided then normalization is not performed.
 */
public class SingleThreadReadConverter<T extends Spot> implements Converter {
    List<ReadReader> readers = new ArrayList<>();
    final List<InputStream> istreams = new ArrayList<>();
    final ReadWriter<Read, T> writer;
    final Long readLimit;

    long readCount = 0, baseCount = 0;

    List<Integer> turnsList = new ArrayList<>();
    int turn;

    public SingleThreadReadConverter(List<InputStream> istreams, ReadWriter<Read, T> writer) {
        this.writer = writer;

        for (int readerIndex = 0; readerIndex < istreams.size(); readerIndex++) {
            this.istreams.add(istreams.get(readerIndex));
            readers.add(new ReadReader(String.valueOf(readerIndex + 1)));
            turnsList.add(readerIndex);
        }

        readLimit = 0L;
    }

    public SingleThreadReadConverter(List<InputStream> istreams, ReadWriter<Read, T> writer, Long readLimit) {
        this.writer = writer;

        for (int readerIndex = 0; readerIndex < istreams.size(); readerIndex++) {
            this.istreams.add(istreams.get(readerIndex));
            readers.add(new ReadReader(String.valueOf(readerIndex + 1)));
            turnsList.add(readerIndex);
        }

        this.readLimit = readLimit;
    }

    public SingleThreadReadConverter(
            List<InputStream> istreams,
            List<QualityNormalizer> normalizers,
            ReadWriter<Read, T> writer,
            Long readLimit) {

        this.writer = writer;

        for (int readerIndex = 0; readerIndex < istreams.size(); readerIndex++) {
            this.istreams.add(istreams.get(readerIndex));
            readers.add(new ReadReader(normalizers.get(normalizers.size() == istreams.size() ? readerIndex : 0), String.valueOf(readerIndex + 1)));
            turnsList.add(readerIndex);
        }

        this.readLimit = readLimit;
    }

    private boolean keepRunning() {
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
                for (int yield = YIELD_CYCLES_FOR_ERROR_CHECKING; yield > 0 && keepRunning(); --yield) {
                    writer.write(convert());
                }

                if (!writer.isOk()) {
                    throw new ConverterPanicException();
                }
            } while (!isDone());
        } catch (ConverterEOFException ignored) {
        } catch (Exception e) {
            if (e instanceof ReadWriterMemoryLimitException) {
                throw e;
            } else {
                throw new RuntimeException(e);
            }
        }
    }

    public void runOnce() {
        try {
            if (!isDone()) {
                writer.write(convert());

                if (!writer.isOk()) {
                    throw new ConverterPanicException();
                }
            }
        } catch (ConverterEOFException ignored) {
        } catch (Exception e) {
            if (e instanceof ReadWriterMemoryLimitException) {
                throw e;
            } else {
                throw new RuntimeException(e);
            }
        }
    }

    public boolean isDone() {
        return (turnsList.isEmpty() || !keepRunning());
    }

    private Read convert() {
        try {
            if (turn >= turnsList.size()) {
                turn = 0;
            }

            Read spot = readers.get(turnsList.get(turn)).read(istreams.get(turnsList.get(turn)));

            turn++;

            readCount++;
            baseCount += spot.getBaseCount();

            return spot;
        } catch(EOFException e){
            turnsList.remove(turn);
            throw new ConverterEOFException(readCount);
        } catch(ConverterException e){
            throw e;
        } catch(Exception cause) {
            throw new ConverterException(cause);
        }
    }
}
