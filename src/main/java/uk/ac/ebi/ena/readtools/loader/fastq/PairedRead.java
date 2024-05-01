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

import java.util.StringJoiner;
import uk.ac.ebi.ena.readtools.loader.common.writer.Spot;

/** Holds unpaired or paired read information. */
public class PairedRead implements Spot {

  public final String name;

  public final Read forward;
  public final Read reverse;

  public PairedRead(String name, Read read) {
    if (name == null || name.isEmpty()) {
      throw new IllegalArgumentException("Invalid read name.");
    }

    if (read == null) {
      throw new IllegalArgumentException("Read cannot be null.");
    }

    this.forward = read;
    this.reverse = null;
    this.name = name;
  }

  /**
   * Constructor for a paired read.
   *
   * @param name
   * @param forward
   * @param reverse
   */
  public PairedRead(String name, Read forward, Read reverse) {
    if (name == null || name.isEmpty()) {
      throw new IllegalArgumentException("Invalid read name.");
    }

    if (forward == null && reverse == null) {
      throw new IllegalArgumentException("Both reads cannot be null.");
    }

    this.forward = forward;
    this.reverse = reverse;

    this.name = name;
  }

  @Override
  public String getName() {
    return name;
  }

  @Override
  public long getBaseCount() {
    long fwd = forward != null ? forward.getBaseCount() : 0;
    long rev = reverse != null ? reverse.getBaseCount() : 0;

    return fwd + rev;
  }

  public boolean isPaired() {
    return reverse != null && forward != null;
  }

  /**
   * Convenient method for getting the non-null unpaired read.
   *
   * @return the non-null unpaired read. Returns null if this is a paired read and both forward and
   *     reverse reads are present.
   */
  public Read getUnpaired() {
    if (forward == null) {
      return reverse;
    } else if (reverse == null) {
      return forward;
    } else {
      return null;
    }
  }

  @Override
  public String toString() {
    return new StringJoiner(", ", PairedRead.class.getSimpleName() + "[", "]")
        .add("name='" + name + "'")
        .add("first=" + forward)
        .add("second=" + reverse)
        .toString();
  }

  @Override
  public long getSizeBytes() {
    return name.getBytes().length + forward.getSizeBytes() + reverse.getSizeBytes();
  }
}
