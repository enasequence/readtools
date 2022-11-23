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
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicReference;

import uk.ac.ebi.ena.readtools.common.reads.QualityNormalizer;
import uk.ac.ebi.ena.readtools.loader.common.FileCompression;
import uk.ac.ebi.ena.readtools.loader.common.converter.ReadConverter;
import uk.ac.ebi.ena.readtools.loader.common.converter.SingleThreadReadConverter;
import uk.ac.ebi.ena.readtools.loader.common.writer.ReadWriter;
import uk.ac.ebi.ena.readtools.loader.common.writer.ReadWriterException;
import uk.ac.ebi.ena.readtools.loader.common.writer.Spot;
import uk.ac.ebi.ena.readtools.loader.fastq.FastqIterativeWriter.READ_TYPE;

public class
FastqIterativeWriterIterator implements Iterator<PairedRead>, ReadWriter<PairedRead, Spot> {
    /** This used to be SynchronousQueue but was replaced due to its high CPU overhead. This one has overhead as well
     * (possibly due to slow consumer), but it is much lower. Increasing capacity helps but not massively. */
    private final BlockingQueue<PairedRead> queue = new ArrayBlockingQueue<>(1_000);
    private final AtomicReference<PairedRead> currentElement = new AtomicReference<>();
    private boolean wasCascadeErrors = false;

    FastqIterativeWriterIterator(File tmp_folder,
                                 int spill_page_size, //only for paired
                                 long spill_page_size_bytes, //only for paired
                                 long spill_abandon_limit_bytes,
                                 READ_TYPE read_type,
                                 File[] files,
                                 final QualityNormalizer normalizers[],
                                 Long readLimit) throws SecurityException, IOException {

        ReadWriter<Read, PairedRead> consumer;

        switch (read_type) {
            case SINGLE:
                consumer = new SingleFastqConsumer();
                break;
            case PAIRED:
                consumer = new PairedFastqWriter(
                        tmp_folder, spill_page_size, spill_page_size_bytes, spill_abandon_limit_bytes);
                break;
            default:
                throw new UnsupportedOperationException();
        }

        consumer.setWriter(this);

        if (files.length == 1) {
            ReadConverter converter;
            if (normalizers != null) {
                converter = new ReadConverter(
                        FileCompression.open(files[0]), readLimit, normalizers[0], "1");
            } else {
                converter = new ReadConverter(FileCompression.open(files[0]), readLimit, "1");
            }

            converter.setWriter(consumer);
            converter.run();
        } else {
            List<InputStream> istreams = new ArrayList<>();
            for (File file : files) {
                istreams.add(FileCompression.open(file));
            }

            SingleThreadReadConverter<PairedRead> converter = (normalizers != null)
                    ? new SingleThreadReadConverter<>(istreams, Arrays.asList(normalizers), consumer, readLimit)
                    : new SingleThreadReadConverter<>(istreams, consumer, readLimit);
            converter.run();
        }
    }

    @Override
    public void
    cascadeErrors() throws ReadWriterException {
        wasCascadeErrors = true;
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
        //There is a flaw in this method. If the producer stops due to an error, while the current thread is in the
        //loop below, then a 'false' value will be returned from this method telling the user of the iterator that
        //there is no element left causing iteration to stop normally. But the right way to end would be to
        //throw the exception that was caught in 'storedException' so the user can know that execution actually failed.
        //TODO update this class so that end of data is detected reliably and errors are thrown and caught properly.
        try {
            boolean currentPairedReadUpdated;
            do {
                //The timeout below is used as way to determine that end of data has been reached.
                //This used to be 1 second but was increased because some unit tests kept randomly failing
                //due to timeout occurring before the data got added to the queue.
                //This is most probably due to converter being slow to start generating the data.
                PairedRead newPairedRead = queue.poll(2L, TimeUnit.SECONDS);

                currentPairedReadUpdated = currentElement.compareAndSet(null, newPairedRead);
            } while (!wasCascadeErrors && !currentPairedReadUpdated);
        } catch (InterruptedException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
        return currentElement.get() != null;
    }

    @Override
    public PairedRead
    next() {
        return currentElement.getAndSet(null);
    }

    @Override
    public void
    remove() {
        throw new RuntimeException("N07 iMPl3m3nt3D");
    }

    @Override
    public boolean
    isOk() {
        return true;
    }
}
