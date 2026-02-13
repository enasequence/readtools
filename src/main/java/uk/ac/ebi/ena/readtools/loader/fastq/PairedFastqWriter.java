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

import static uk.ac.ebi.ena.readtools.loader.common.writer.ReadWriterException.ErrorType.INVALID_READ_NAME;

import java.io.File;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import uk.ac.ebi.ena.readtools.common.reads.CasavaRead;
import uk.ac.ebi.ena.readtools.loader.common.writer.AbstractPagedReadWriter;
import uk.ac.ebi.ena.readtools.loader.common.writer.ReadWriterException;

public class PairedFastqWriter extends AbstractPagedReadWriter<Read, PairedRead> {
  // Provided readname structure is @{readkey}{separator:1(.|/|:|_)}{index:1(0:1:2)}
  static final Pattern SPLIT_REGEXP = Pattern.compile("^(.*)(?:[\\.|:|/|_])([0-9]+)$");
  /*
  @ Each sequence identifier line starts with @
  1    <instrument> Characters
  allowed:
  a-z, A-Z, 0-9 and
  underscore
  2    Instrument ID
  <run number> Numerical Run number on instrument
  3    <flowcell
  ID>
  Characters
  allowed:
  a-z, A-Z, 0-9
  4    <lane> Numerical Lane number
  5    <tile> Numerical Tile number
  6    <x_pos> Numerical X coordinate of cluster
  7    <y_pos> Numerical Y coordinate of cluster
  */
  //    A00953:544:HMTFHDSX3:2:1101:6768:1
  //             1        :  2   :    3       :   4  :  5   :   6   :  7
  //    "^([a-zA-Z0-9_-]+:[0-9]+:[a-zA-Z0-9]+:[0-9]+:[0-9]+:[0-9-]+:[0-9-]+)$"
  static final Pattern CASAVA_LIKE_EXCLUDE_REGEXP =
      Pattern.compile("^([a-zA-Z0-9_-]+:[0-9]+:[a-zA-Z0-9_-]+:[0-9]+:[0-9]+:[0-9-]+:[0-9-]+)$");

  public static final int KEY = 1;
  public static final int INDEX = 2;

  Integer index1 = null, index2 = null;

  public PairedFastqWriter(
      File tmp_root,
      int spill_page_size,
      long spill_page_size_bytes,
      long spill_abandon_limit_bytes) {
    super(tmp_root, spill_page_size, spill_page_size_bytes, spill_abandon_limit_bytes);
  }

  public static String getReadKey(String readname) throws ReadWriterException {
    return getReadPart(readname, KEY);
  }

  public static String getPairNumber(String readname) throws ReadWriterException {
    return getReadPart(readname, INDEX);
  }

  private static String getReadPart(String readname, int group) throws ReadWriterException {
    // Try Casava 1.8 format first
    String casavaResult =
        group == KEY
            ? CasavaRead.getBaseNameOrNull(readname)
            : CasavaRead.getReadIndexOrNull(readname);
    if (casavaResult != null) return casavaResult;

    // Fall back to generic split logic for non-Casava names
    Matcher casavaLikeMatcher = CASAVA_LIKE_EXCLUDE_REGEXP.matcher(readname);
    if (!casavaLikeMatcher.find()) {
      Matcher m = SPLIT_REGEXP.matcher(readname);
      if (m.find()) {
        return m.group(group);
      }
    }
    throw new ReadWriterException(
        String.format("Readname [%s] does not match regexp", readname),
        ReadWriterException.ErrorType.INVALID_READ_NAME);
  }

  @Override
  public String getKey(Read spot) {
    try {
      return getReadKey(spot.name);
    } catch (ReadWriterException e) {
      if (INVALID_READ_NAME.equals(e.getErrorType())) {
        return spot.name;
      } else {
        throw e;
      }
    }
  }

  public List<Read> newListBucket() {
    List<Read> list = super.newListBucket();
    list.add(null);
    list.add(null);
    return list;
  }

  public void append(List<Read> list, Read spot) throws ReadWriterException {
    int readIndex = getReadIndex(spot);

    // Initialize index1 and index2 if they are not set yet
    if (index1 == null) {
      index1 = readIndex;
    } else if (index2 == null && readIndex != index1) {
      index2 = readIndex;
    }

    // Check if the readIndex is valid
    if (readIndex != index1 && readIndex != index2) {
      throw new ReadWriterException(
          "Unexpected read pair number: "
              + readIndex
              + "; pair numbers "
              + index1
              + " and "
              + index2
              + " were found previously in the file.",
          ReadWriterException.ErrorType.UNEXPECTED_PAIR_NUMBER);
    }

    int mappedIndex = (readIndex == index1) ? 0 : 1;

    if (list.get(mappedIndex) == null) {
      list.set(mappedIndex, spot);
    } else {
      throw new ReadWriterException(
          "Got same spot twice: " + spot, ReadWriterException.ErrorType.SPOT_DUPLICATE);
    }

    // Check if the list contains any nulls
    if (!list.contains(null)) {
      // Verify and sort the list based on readIndex
      try {
        list.sort(
            (read1, read2) -> {
              int readIndex1 = getReadIndex(read1);
              int readIndex2 = getReadIndex(read2);
              return Integer.compare(readIndex1, readIndex2);
            });
      } catch (RuntimeException e) {
        throw new ReadWriterException(
            "Error sorting reads by read index:" + e.getMessage(),
            ReadWriterException.ErrorType.SORTING_ERROR);
      }
    }
  }

  private int getReadIndex(Read spot) throws ReadWriterException {
    String readIndexStr;
    try {
      readIndexStr = getPairNumber(spot.name);
    } catch (ReadWriterException e) {
      if (INVALID_READ_NAME.equals(e.getErrorType())) {
        readIndexStr = spot.getDefaultReadIndex();
      } else {
        throw e;
      }
    }
    return Integer.parseInt(readIndexStr) - 1;
  }

  @Override
  public PairedRead assemble(final String key, List<Read> list) throws ReadWriterException {
    PairedRead spot =
        list.size() == 1
            ? new PairedRead(key, list.get(0))
            : new PairedRead(key, list.get(0), list.get(1));

    return spot;
  }

  @Override
  public boolean isCollected(List<Read> list) {
    return null != list.get(0) && null != list.get(1);
  }

  @Override
  public PairedRead handleErrors(final String key, List<Read> list) throws ReadWriterException {
    return assemble(key, list);
  }
}
