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
package uk.ac.ebi.ena.readtools.webin.cli.rawreads;

import java.util.Set;
import java.util.concurrent.atomic.AtomicLong;

import uk.ac.ebi.ena.readtools.loader.common.writer.ReadWriter;
import uk.ac.ebi.ena.readtools.loader.common.writer.ReadWriterException;
import uk.ac.ebi.ena.readtools.loader.common.writer.Spot;
import uk.ac.ebi.ena.readtools.loader.fastq.PairedFastqWriter;
import uk.ac.ebi.ena.readtools.loader.fastq.Read;

public class FastqReadScanner implements ReadWriter<Read, Spot> {
    private final String streamName;
    private final Set<String> labels;
    private final AtomicLong count;
    private final BloomWrapper pairingBloomWrapper;
    private final BloomWrapper duplicationsBloomWrapper;
    private final FastqScanner fastqScanner;
    private final int maxLabelSetSize;
    private final int printFreq;

    public FastqReadScanner(
            String streamName, Set<String> labels, AtomicLong count,
            BloomWrapper pairingBloomWrapper, BloomWrapper duplicationsBloomWrapper,
            FastqScanner fastqScanner, int maxLabelSetSize, int printFreq)
    {
        this.streamName = streamName;
        this.labels = labels;
        this.count = count;
        this.pairingBloomWrapper = pairingBloomWrapper;
        this.duplicationsBloomWrapper = duplicationsBloomWrapper;
        this.fastqScanner = fastqScanner;
        this.maxLabelSetSize = maxLabelSetSize;
        this.printFreq = printFreq;
    }

    @Override
    public void cascadeErrors() throws ReadWriterException {
    }

    @Override
    public void write(Read spot) throws ReadWriterException {
        String readNameWithoutPairNumber;
        String pairNumber;

        try {
            readNameWithoutPairNumber = PairedFastqWriter.getReadKey(spot.name);
            pairNumber = PairedFastqWriter.getPairNumber(spot.name);
        } catch (ReadWriterException ignored) {
            readNameWithoutPairNumber = spot.name;
            pairNumber = streamName;
        }

        if (labels.size() < maxLabelSetSize) {
            labels.add(pairNumber);
        }

        count.incrementAndGet();
        pairingBloomWrapper.add(readNameWithoutPairNumber);
        duplicationsBloomWrapper.add(spot.name);

        if (0 == count.get() % printFreq)
            fastqScanner.logProcessedReadNumber(count.get());
    }

    @Override
    public void setWriter(ReadWriter readWriter) {
        throw new RuntimeException("Not implemented");
    }
}
