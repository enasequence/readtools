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
import java.io.InputStream;
import java.util.function.Supplier;

import uk.ac.ebi.ena.readtools.loader.common.writer.ReadWriter;
import uk.ac.ebi.ena.readtools.loader.common.writer.Spot;

public abstract class
AbstractReadConverter<T extends Spot> extends Thread implements Converter<T> {
    private static final int YIELD_CYCLES_FOR_ERROR_CHECKING = 362;//16384;

    protected final InputStream  istream;

    protected final Long readLimit;

    private volatile long readCount = 0, baseCount = 0;

    protected volatile ReadWriter<T, ?> readWriter;
    protected volatile boolean      is_ok = true;
    protected volatile Throwable    stored_exception;
    
    protected AbstractReadConverter(InputStream istream) {
        this(istream, null);
    }

    /**
     *
     * @param istream
     * @param readLimit Only read limited amount of reads.
     */
    protected AbstractReadConverter(InputStream istream, Long readLimit) {
        this.istream = new BufferedInputStream( istream, 1024 * 1024 );
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
    
    public void setWriter(ReadWriter<T, ?> writer) {
        this.readWriter = writer;
    }
    
    public boolean isOk() {
        return is_ok;
    }
    
    public Throwable getStoredException() {
        return stored_exception;
    }

    //TODO: scheduler should be fair and based on different principle
    public final void run() {
        try {
            Supplier<Boolean> keepRunning = createKeepRunning();

            begin();

            do {
                synchronized (readWriter) {
                    for (int yield = YIELD_CYCLES_FOR_ERROR_CHECKING; yield > 0 && keepRunning.get(); --yield)
                        readWriter.write(convert());
                }

                if (!readWriter.isOk())
                    throw new ConverterPanicException();

            } while (keepRunning.get());
        } catch (ConverterEOFException ignored) {
        } catch (ConverterPanicException e) {
            is_ok = false;
            this.stored_exception = e;
        } catch (Throwable e) {
            this.stored_exception = e;
            is_ok = false;
        } finally {
            end();
        }
    }

    protected void begin() {}

    protected void end() {}

    private Supplier<Boolean> createKeepRunning() {
        if (readLimit == null) {
            return () -> true;
        } else {
            return () -> readCount < readLimit;
        }
    }

    //Re-implement if you need special type of feeding
    private T convert() {
        T spot = null;

        try {
            spot = convert(istream);
            ++readCount;
            baseCount += spot.getBaseCount();


            return spot;
        } catch( EOFException e ){
            throw new ConverterEOFException(readCount);
        } catch( ConverterException e ){
            throw e;
        } catch( Throwable cause ) {
            throw new ConverterException(cause);
        }
    }
}
