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
import java.util.*;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicReference;

import uk.ac.ebi.ena.readtools.common.reads.QualityNormalizer;
import uk.ac.ebi.ena.readtools.loader.common.FileCompression;
import uk.ac.ebi.ena.readtools.loader.common.converter.Converter;
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
    private final Converter converter;
    private final Queue<PairedRead> queue = new LinkedList<>();

    FastqIterativeWriterIterator(File tmp_folder,
                                 int spill_page_size, //only for paired
                                 long spill_page_size_bytes, //only for paired
                                 long spill_abandon_limit_bytes,
                                 READ_TYPE read_type,
                                 File[] files,
                                 final QualityNormalizer[] normalizers,
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
            if (normalizers != null) {
                converter = new ReadConverter(
                        FileCompression.open(files[0]), consumer, readLimit, normalizers[0], "1");
            } else {
                converter = new ReadConverter(FileCompression.open(files[0]), consumer, readLimit, "1");
            }
        } else {
            List<InputStream> istreams = new ArrayList<>();
            for (File file : files) {
                istreams.add(FileCompression.open(file));
            }

            converter = (normalizers != null)
                    ? new SingleThreadReadConverter<>(istreams, Arrays.asList(normalizers), consumer, readLimit)
                    : new SingleThreadReadConverter<>(istreams, consumer, readLimit);
        }
    }

    @Override
    public void
    cascadeErrors() {
    }

    @Override
    public void
    write(PairedRead spot) throws ReadWriterException {
        queue.add(spot);
    }

    @Override
    public void
    setWriter(ReadWriter<Spot, ?> readWriter) {
        throw new RuntimeException("N07 iMPl3m3nt3D");
    }

    @Override
    public boolean
    hasNext() {
        iterate();
        return !queue.isEmpty();
    }

    @Override
    public PairedRead
    next() {
        iterate();
        return queue.poll();
    }

    private void
    iterate() {
        while (!converter.isDone() && queue.isEmpty()) {
            converter.runOnce();
        }
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
