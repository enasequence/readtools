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
package uk.ac.ebi.ena.readtools.loader.fastq;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicReference;

import uk.ac.ebi.ena.readtools.common.reads.QualityNormalizer;
import uk.ac.ebi.ena.readtools.loader.common.FileCompression;
import uk.ac.ebi.ena.readtools.loader.common.converter.ConverterException;
import uk.ac.ebi.ena.readtools.loader.common.converter.ReadConverter;
import uk.ac.ebi.ena.readtools.loader.common.writer.ReadWriter;
import uk.ac.ebi.ena.readtools.loader.common.writer.ReadWriterException;
import uk.ac.ebi.ena.readtools.loader.common.writer.Spot;
import uk.ac.ebi.ena.readtools.loader.fastq.FastqIterativeWriter.READ_TYPE;

public class
FastqIterativeWriterIterator implements Iterator<PairedRead>, ReadWriter<PairedRead, Spot> {

    private BlockingQueue<PairedRead> queue = new SynchronousQueue<PairedRead>();
    private AtomicReference<PairedRead> current_element = new AtomicReference<PairedRead>();
    private AtomicBoolean was_cascade_errors = new AtomicBoolean(false);
    private Exception storedException;

    FastqIterativeWriterIterator(File tmp_folder,
                                 int spill_page_size, //only for paired
                                 long spill_page_size_bytes, //only for paired
                                 long spill_abandon_limit_bytes,
                                 READ_TYPE read_type,
                                 File[] files,
                                 final QualityNormalizer normalizers[]) throws SecurityException, IOException {
        ReadWriter<Read, PairedRead> writer = null;

        switch (read_type) {
            case SINGLE:
                writer = new SingleFastqConsumer();
                break;
            case PAIRED:
                writer = new PairedFastqWriter(tmp_folder, spill_page_size, spill_page_size_bytes, spill_abandon_limit_bytes);
                break;
            default:
                throw new UnsupportedOperationException();
        }

        writer.setWriter(this);

        ArrayList<ReadConverter> converters = new ArrayList<>();

        ExecutorService executorService = Executors.newCachedThreadPool();
        List<Future<?>> futures = Collections.synchronizedList(new ArrayList<>());

        int attr = 1;
        for (File file : files) {
            final String default_attr = Integer.toString(attr++);
            final int nindex = normalizers.length == files.length ? attr - 2 : 0;

            ReadConverter converter = new ReadConverter(FileCompression.open(file), normalizers[nindex], default_attr);
            converter.setWriter(writer);
            converters.add(converter);

            futures.add(executorService.submit(converter::run));
        }

        ReadWriter<Read, PairedRead> finalWriter = writer;
        executorService.submit(() -> lifecycle(converters, futures, finalWriter));
    }

    private void
    lifecycle(final ArrayList<ReadConverter> converters,
              final List<Future<?>> futures,
              final ReadWriter<?, ?> writer) {

        try {
            futures.forEach(f -> {
                try {
                    f.get();
                } catch (InterruptedException ignored) {
                } catch (ExecutionException e) {
                    throw new ConverterException(e);
                }
            });

            for (ReadConverter converter : converters) {
                if (!converter.isOk()) {
                    throw new ConverterException(converter.getStoredException());
                }
            }
            writer.cascadeErrors();

        } catch (Exception dfe) {
            storedException = dfe;
        }
    }

    @Override
    public void
    cascadeErrors() throws ReadWriterException {
        was_cascade_errors.lazySet(true);
    }

    @Override
    public void
    write(PairedRead spot) throws ReadWriterException {
        try {
            queue.put(spot);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
    }

    @Override
    public void
    setWriter(ReadWriter<Spot, ?> readWriter) {
        throw new RuntimeException("N07 iMPl3m3nt3D");
    }

    @Override
    public boolean
    hasNext() {
        if (null != storedException) {
            throw new RuntimeException(storedException);
        }

        try {
            while (!was_cascade_errors.get()
                    && false == current_element
                    .compareAndSet(null, queue.poll(1L, TimeUnit.SECONDS))) {
                ;
            }
        } catch (InterruptedException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
        return null != current_element.get();
    }

    @Override
    public PairedRead
    next() {
        if (null != storedException) {
            throw new RuntimeException(storedException);
        }

        return current_element.getAndSet(null);
    }

    @Override
    public void
    remove() {
        throw new RuntimeException("N07 iMPl3m3nt3D");
    }

    public Throwable
    getStoredException() {
        return storedException;
    }

    @Override
    public boolean
    isOk() {
        return true;
    }
}
