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
import uk.ac.ebi.ena.readtools.loader.common.writer.ReadWriterException;
import uk.ac.ebi.ena.readtools.loader.fastq.PairedFastqWriter;
import uk.ac.ebi.ena.readtools.loader.fastq.Read;

public class FastqReadScanner extends InsdcStandardCheckingScanner {
  private final String streamName;
  private final Set<String> labels;
  private final BloomWrapper pairingBloomWrapper;
  private final BloomWrapper duplicationsBloomWrapper;
  private final int maxLabelSetSize;

  public FastqReadScanner(
      String streamName,
      Set<String> labels,
      BloomWrapper pairingBloomWrapper,
      BloomWrapper duplicationsBloomWrapper,
      int maxLabelSetSize,
      int printFreq) {
    super(printFreq);

    this.streamName = streamName;
    this.labels = labels;
    this.pairingBloomWrapper = pairingBloomWrapper;
    this.duplicationsBloomWrapper = duplicationsBloomWrapper;
    this.maxLabelSetSize = maxLabelSetSize;
  }

  @Override
  protected void logProcessedReadNumber(long cnt) {
    ;
  }

  public void write(Read read) throws ReadWriterException {
    super.write(read);

    String readNameWithoutPairNumber;
    String pairNumber;

    try {
      readNameWithoutPairNumber = PairedFastqWriter.getReadKey(read.getName());
      pairNumber = PairedFastqWriter.getPairNumber(read.getName());
    } catch (ReadWriterException ignored) {
      readNameWithoutPairNumber = read.getName();
      pairNumber = streamName;
    }

    if (labels.size() < maxLabelSetSize) {
      labels.add(pairNumber);
    }

    pairingBloomWrapper.add(readNameWithoutPairNumber);
    duplicationsBloomWrapper.add(read.getName());
  }
}
