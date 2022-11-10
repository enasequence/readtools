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

import uk.ac.ebi.ena.readtools.loader.fastq.PairedFastqWriter;
import uk.ac.ebi.ena.readtools.loader.fastq.Read;

import java.io.EOFException;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.List;

/**
 * Similar to {@link AutoNormalizeQualityReadConverter}, but here, base quality normalizer is provided explicitly. If no normalizer
 * is provided then normalization is not performed.
 */
public class SingleThreadReadConverter {
    static final int YIELD_CYCLES_FOR_ERROR_CHECKING = 362;
    List<ReadReader> readers = new ArrayList<>();
    List<InputStream> istreams = new ArrayList<>();
    PairedFastqWriter writer;

    long readCount = 0, baseCount = 0;

    List<Integer> turnsList = new ArrayList<>();
    int turn;
    private Exception storedException;

    public SingleThreadReadConverter(InputStream... istreams) {
        for (int readerIndex = 0; readerIndex < istreams.length; readerIndex++) {
            this.istreams.add(istreams[readerIndex]);
            readers.add(new ReadReader(String.valueOf(readerIndex + 1)));
            turnsList.add(readerIndex);
        }
    }

    public long getReadCount() {
        return readCount;
    }

    public long getBaseCount() {
        return baseCount;
    }

    public Exception getStoredException() {
        return storedException;
    }

    public final void run() {
        try {
            do {
                for (int yield = YIELD_CYCLES_FOR_ERROR_CHECKING; yield > 0; --yield) {
                    writer.write(convert());
                }

                if (!writer.isOk()) {
                    throw new ConverterPanicException();
                }
            } while (turnsList.size() > 0);
        } catch (ConverterEOFException ignored) {
        } catch (Exception e) {
            this.storedException = e;
        }
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
