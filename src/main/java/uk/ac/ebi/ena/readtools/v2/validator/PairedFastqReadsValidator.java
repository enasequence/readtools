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
package uk.ac.ebi.ena.readtools.v2.validator;

import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import uk.ac.ebi.ena.readtools.v2.read.FastqRead;
import uk.ac.ebi.ena.readtools.webin.cli.rawreads.BloomWrapper;

public class PairedFastqReadsValidator extends FastqReadsValidator {
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
  private static final Pattern CASAVA_LIKE_EXCLUDE_REGEXP =
      Pattern.compile("^([a-zA-Z0-9_-]+:[0-9]+:[a-zA-Z0-9_-]+:[0-9]+:[0-9]+:[0-9-]+:[0-9-]+)$");
  // Provided readname structure is @{readkey}{separator:1(.|/|:|_)}{index:1(0:1:2)}
  private static final Pattern SPLIT_REGEXP = Pattern.compile("^(.*)(?:[\\.|:|/|_])([1234])$");

  private final String providerName;
  private final BloomWrapper pairingBloomWrapper;
  private final Set<String> labels;

  private long addCount = 0;

  public PairedFastqReadsValidator(
      long readCountLimit,
      String providerName,
      BloomWrapper pairingBloomWrapper,
      Set<String> labels) {
    super(readCountLimit);
    this.providerName = providerName;
    this.pairingBloomWrapper = pairingBloomWrapper;
    this.labels = labels;
  }

  public long getAddCount() {
    return addCount;
  }

  @Override
  protected void extraReadsValidation(ReadStyle readStyle, long readCount, FastqRead read)
      throws ReadsValidationException {
    addCount++;

    if (pairingBloomWrapper != null && labels != null) {
      if (readStyle == ReadStyle.CASAVA18) {
        pairingBloomWrapper.add(getCasavaReadNameWithoutIndex(read.getName(), readCount));
        labels.add(getCasavaReadIndex(read.getName(), readCount));
      } else {
        pairingBloomWrapper.add(getNonCasavaReadNameWithoutIndex(read.getName()));
        labels.add(getNonCasavaReadIndex(read.getName()));
      }
    }
  }

  private String getCasavaReadNameWithoutIndex(String readName, long readIndex)
      throws ReadsValidationException {
    Matcher matcher = P_CASAVA_18_NAME.matcher(readName);
    if (!matcher.matches()) {
      throw new ReadsValidationException(
          String.format("Line [%s] does not match %s regexp", readName, ReadStyle.CASAVA18),
          readIndex);
    }
    return matcher.group(1);
  }

  private String getCasavaReadIndex(String readName, long readIndex)
      throws ReadsValidationException {
    Matcher matcher = P_CASAVA_18_NAME.matcher(readName);
    if (!matcher.matches()) {
      throw new ReadsValidationException(
          String.format("Line [%s] does not match %s regexp", readName, ReadStyle.CASAVA18),
          readIndex);
    }
    return matcher.group(3);
  }

  private String getNonCasavaReadNameWithoutIndex(String readName) {
    Matcher casavaLikeMatcher = CASAVA_LIKE_EXCLUDE_REGEXP.matcher(readName);
    if (!casavaLikeMatcher.find()) {
      Matcher m = SPLIT_REGEXP.matcher(readName);
      if (m.find()) {
        return m.group(1);
      }
    }
    return readName;
  }

  private String getNonCasavaReadIndex(String readName) {
    Matcher casavaLikeMatcher = CASAVA_LIKE_EXCLUDE_REGEXP.matcher(readName);
    if (!casavaLikeMatcher.find()) {
      Matcher m = SPLIT_REGEXP.matcher(readName);
      if (m.find()) {
        return m.group(2);
      }
    }
    return providerName;
  }
}
