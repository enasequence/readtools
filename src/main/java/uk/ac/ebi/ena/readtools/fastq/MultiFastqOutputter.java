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
package uk.ac.ebi.ena.readtools.fastq;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.Log;
import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

/**
 * Note<br>
 * Makes use of a collection based caching implementation that may produce results in which a paired
 * read is incorrectly written to the unpaired file if the two mates are located far away from each
 * other in the collection.
 */
public class MultiFastqOutputter {
  private static final Log log = Log.getInstance(MultiFastqOutputter.class);
  private Map<FastqRead, FastqRead> readSet = new TreeMap<FastqRead, FastqRead>();
  private int maxCacheSize =
      Integer.parseInt(System.getProperty("fastq-dumper.cache-size", Integer.toString(100000)));

  private static Comparator<FastqRead> byGenerationComparator =
      (o1, o2) -> (int) (o1.generation - o2.generation);

  private long generation = 0;
  private OutputStream[] streams;

  private SAMFileHeader headerForOverflowWriter;
  private SAMFileWriter writer;
  private OutputStream cacheOverFlowStream;
  private byte[] prefix;
  private long counter = 1;
  private SAMFileHeader header;

  List<FastqRead> list = new ArrayList<FastqRead>();

  public MultiFastqOutputter(
      OutputStream[] streams, OutputStream cacheOverFlowStream, SAMFileHeader header) {
    this.streams = streams;
    this.cacheOverFlowStream = cacheOverFlowStream;
    this.header = header;
  }

  public byte[] getPrefix() {
    return prefix;
  }

  public void setPrefix(byte[] prefix) {
    this.prefix = prefix;
  }

  public long getCounter() {
    return counter;
  }

  protected void write(FastqRead read, OutputStream stream) throws IOException {
    if (prefix == null) {
      stream.write(read.data);
    } else {
      streams[read.templateIndex].write('@');
      streams[read.templateIndex].write(prefix);
      streams[read.templateIndex].write('.');
      streams[read.templateIndex].write(String.valueOf(counter).getBytes());
      streams[read.templateIndex].write(' ');
      streams[read.templateIndex].write(read.data, 1, read.data.length - 1);
    }
  }

  protected void foundCollision(FastqRead read) {
    FastqRead anchor = readSet.remove(read);
    try {
      write(anchor, streams[anchor.templateIndex]);
      write(read, streams[read.templateIndex]);
      counter++;
    } catch (IOException e) {
      throw new RuntimeException(e);
    }
  }

  protected void kickedFromCache(FastqRead read) {
    if (writer == null) {
      log.info("Creating overflow BAM file.");
      headerForOverflowWriter = header.clone();
      headerForOverflowWriter.setSortOrder(SAMFileHeader.SortOrder.queryname);

      writer =
          new SAMFileWriterFactory()
              .makeBAMWriter(headerForOverflowWriter, false, cacheOverFlowStream);
    }
    SAMRecord r = read.toSAMRecord(writer.getFileHeader());
    writer.addAlignment(r);
  }

  protected void purgeCache() {
    long time1 = System.nanoTime();
    list.clear();
    for (FastqRead read : readSet.keySet()) list.add(read);

    Collections.sort(list, byGenerationComparator);
    for (int i = 0; i < list.size() / 2; i++) {
      readSet.remove(list.get(i));
      kickedFromCache(list.get(i));
    }

    list.clear();
    long time2 = System.nanoTime();
    log.debug(String.format("Cache purged in %.2fms.\n", (time2 - time1) / 1000000f));
  }

  public void writeRead(byte[] name, int flags, byte[] bases, byte[] scores) {
    FastqRead read =
        new FastqRead(bases.length, name, true, getSegmentIndexInTemplate(flags), bases, scores);
    read.generation = generation++;
    if (read.templateIndex == 0) {
      try {
        write(read, streams[0]);
        counter++;
      } catch (IOException e) {
        throw new RuntimeException(e);
      }
      return;
    }

    if (readSet.containsKey(read)) {
      foundCollision(read);
    } else {
      readSet.put(read, read);

      if (readSet.size() > maxCacheSize) purgeCache();
    }
  }

  public void finish() {
    for (FastqRead read : readSet.keySet()) kickedFromCache(read);

    readSet.clear();
    if (writer != null) writer.close();
    writer = null;
  }

  /**
   * For now this is to identify the right buffer to use.
   *
   * @param flags read bit flags
   * @return 0 for non-paired or other rubbish which could not be reliably paired, 1 for first in
   *     pair and 2 for second in pair
   */
  private int getSegmentIndexInTemplate(int flags) {
    if ((flags & 1) == 0) return 0;

    if ((flags & 64) != 0) return 1;
    else return 2;
  }
}
