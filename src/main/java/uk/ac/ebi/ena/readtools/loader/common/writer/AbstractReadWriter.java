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
package uk.ac.ebi.ena.readtools.loader.common.writer;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;


public abstract class
AbstractReadWriter<T1 extends Spot, T2 extends Spot> implements ReadWriter<T1, T2> {
    protected Map<String, List<T1>> spots = null;
    protected ReadWriter<T2, ?> readWriter;

    private long log_time = System.currentTimeMillis();
    private long log_interval = 60 * 1000;
    private long assembled = 0;
    private long ate = 0;
    protected long spotsSizeBytes = 0;
    protected boolean verbose = false;

    public AbstractReadWriter<T1, T2>
    setVerbose(boolean verbose) {
        this.verbose = verbose;
        return this;
    }

    public AbstractReadWriter() {
        this(1024 * 1024);
    }

    protected AbstractReadWriter(int map_size) {
        spots = new HashMap<>(map_size);
    }

    @Override
    public void
    setWriter(ReadWriter<T2, ? extends Spot> readWriter) {
        if (!readWriter.equals(this))
            this.readWriter = readWriter;
    }

    public abstract String
    getKey(T1 spot) throws ReadWriterException;

    public abstract T2
    assemble(final String key, List<T1> list) throws ReadWriterException;

    public void
    append(List<T1> list, T1 spot) throws ReadWriterException {
        list.add(spot);
    }

    public List<T1>
    newListBucket() {
        return new ArrayList<T1>();
    }

    public abstract boolean
    isCollected(List<T1> list);

    public void
    cascadeErrors() throws ReadWriterException {
        for (Entry<String, List<T1>> entry : spots.entrySet()) {
            if (null != readWriter)
                readWriter.write(handleErrors(entry.getKey(), entry.getValue()));
            else
                System.out.println("<?> " + handleErrors(entry.getKey(), entry.getValue()));
        }

        if (null != readWriter)
            readWriter.cascadeErrors();
    }

    public abstract T2
    handleErrors(final String key, List<T1> list) throws ReadWriterException;

    public void
    write(T1 spot) throws ReadWriterException {
        String key = getKey(spot);
        List<T1> bucket;

        if (!spots.containsKey(key)) {
            bucket = newListBucket();
            spots.put(key, bucket);
        } else {
            bucket = spots.get(key);
        }

        ate++;
        append(bucket, spot);
        spotsSizeBytes += spot.getSizeBytes();

        if (isCollected(bucket)) {
            T2 assembly = assemble(key, bucket);
            assembled++;

            List<T1> removed = spots.remove(key);
            spotsSizeBytes -= bucketSize(removed);

            if (null != readWriter) {
                readWriter.write(assembly);
            } else {
                System.out.println(assembly);
            }
        }

        if (verbose) {
            long time = System.currentTimeMillis();

            if (time > log_time) {
                log_time = time + log_interval;
                System.out.println(String.format(
                        "Ate: %d,\tAssembled: %d,\tDiff: %d", ate, assembled, ate - (assembled << 1)));
            }
        }
    }

    private long bucketSize(List<T1> bucket) {
        return bucket.stream().mapToLong(Spot::getSizeBytes).sum();
    }
}
