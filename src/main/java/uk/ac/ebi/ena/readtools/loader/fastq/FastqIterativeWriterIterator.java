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

import uk.ac.ebi.ena.readtools.common.reads.QualityNormalizer;
import uk.ac.ebi.ena.readtools.loader.common.FileCompression;
import uk.ac.ebi.ena.readtools.loader.common.converter.ConverterException;
import uk.ac.ebi.ena.readtools.loader.common.converter.ReadConverter;
import uk.ac.ebi.ena.readtools.loader.common.writer.ReadWriter;
import uk.ac.ebi.ena.readtools.loader.common.writer.ReadWriterException;
import uk.ac.ebi.ena.readtools.loader.common.writer.Spot;
import uk.ac.ebi.ena.readtools.loader.fastq.FastqIterativeWriter.READ_TYPE;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicReference;

public class
FastqIterativeWriterIterator implements Iterator<PairedRead>, ReadWriter<PairedRead, Spot> {
    private static final long CYCLE_TIMEFRAME = 0L;

    /** This used to be SynchronousQueue but was replaced due to its high CPU overhead. This one has overhead as well
     * (possibly due to slow consumer), but it is much lower. Increasing capacity helps but not massively. */
    private final BlockingQueue<PairedRead> queue = new ArrayBlockingQueue<>(1_000);
    private final AtomicReference<PairedRead> current_element = new AtomicReference<PairedRead>();
    private final AtomicBoolean was_cascade_errors = new AtomicBoolean(false);
    private volatile Exception storedException;

    FastqIterativeWriterIterator(File tmp_folder,
                                 int spill_page_size, //only for paired
                                 long spill_page_size_bytes, //only for paired
                                 long spill_abandon_limit_bytes,
                                 READ_TYPE read_type,
                                 File[] files,
                                 final QualityNormalizer normalizers[],
                                 Long readLimit) throws SecurityException, IOException {
        ReadWriter<Read, PairedRead> consumer = null;

        switch (read_type) {
            case SINGLE:
                consumer = new SingleFastqConsumer();
                break;
            case PAIRED:
                consumer = new PairedFastqWriter(tmp_folder, spill_page_size, spill_page_size_bytes, spill_abandon_limit_bytes);
                break;
            default:
                throw new UnsupportedOperationException();
        }

        consumer.setWriter(this);

        ArrayList<ReadConverter> producers = new ArrayList<>();

        int attr = 1;

        for (File file : files) {
            final String default_attr = Integer.toString(attr++);

            //A producer extends Thread class and runs asynchronously.
            //Producer reads data from files and send it here where it gets collected inside a queue.
            //The data in the queue is them read back as the iterator gets used.
            //The producing and the consuming, therefore, happens asynchronously.
            ReadConverter producer;
            if (normalizers != null) {
                int nindex = normalizers.length == files.length ? attr - 2 : 0;

                producer = new ReadConverter(FileCompression.open(file), readLimit, normalizers[nindex], default_attr);
            } else {
                producer = new ReadConverter(FileCompression.open(file), readLimit, default_attr);
            }

            producer.setWriter(consumer);
            producer.setName(file.getPath());
            producer.start();

            producers.add(producer);
        }

        //Start a new thread that monitors the producers created above and collects errors they complete with.
        new Thread(lifecycle(producers, consumer), "lifecycle").start();
    }

    private Runnable
    lifecycle(final ArrayList<ReadConverter> producers,
              final ReadWriter<?, ?> consumer_root) {
        return new Runnable() {
            public void
            run() {
                try {
                    boolean again;
                    do {
                        again = false;
                        for (ReadConverter producer : producers) {
                            if (producer.isAlive()) {
                                try {
                                    producer.join(CYCLE_TIMEFRAME);
                                } catch (InterruptedException ie) {
                                    again = true;
                                    System.out.printf("%s was interrupted\n", producer.getName());
                                }
                            } else if (!producer.isOk()) {
                                throw new ConverterException(producer.getStoredException());
                            }
                        }
                    } while (again);

                    for (ReadConverter producer : producers) {
                        if (!producer.isOk()) {
                            throw new ConverterException(producer.getStoredException());
                        }
                    }

                    consumer_root.cascadeErrors();

                } catch (Exception dfe) {
                    storedException = dfe;
                }
            }
        };
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
            boolean errorsOccurred;
            boolean pairedReadUpdated;
            do {
                errorsOccurred = was_cascade_errors.get();

                //The timeout below is used as way to determine that end of data has been reached.
                //This used to be 1 second but was increased because some unit tests kept randomly failing
                //due to timeout occurring before the data got added to the queue. This is most probably due to converter
                //being slow to start generating the data.
                PairedRead newPairedRead = queue.poll(2L, TimeUnit.SECONDS);

                pairedReadUpdated = current_element.compareAndSet(null, newPairedRead);

            } while (!errorsOccurred && !pairedReadUpdated);
        } catch (InterruptedException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
        return current_element.get() != null;
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
