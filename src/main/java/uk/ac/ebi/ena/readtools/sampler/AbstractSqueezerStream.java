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
package uk.ac.ebi.ena.readtools.sampler;

import java.io.IOException;
import java.io.OutputStream;
import java.util.BitSet;

public abstract class AbstractSqueezerStream {
  OutputStream ostream;
  BitSet remains;

  public void setOutputStream(OutputStream ostream) {
    this.ostream = ostream;
  }

  BitSet writeSet(BitSet set) throws IOException {
    BitSet remains = new BitSet(0);
    int mod = set.size() % 8;
    if (0 != mod) {
      remains = set.get(set.size() - mod, set.size());
    }

    int upto = set.size() - mod;
    int value = 0;
    for (int i = 0; i < upto; ++i) {
      if (set.get(i)) value = value | 1 << (i % 8);

      if (0 == i % 8) {
        ostream.write(value);
        value ^= value;
      }
    }
    return remains;
  }

  public abstract BitSet squeezeMethod(BitSet remains, byte bytes[]);

  public void squeeze(byte bytes[]) throws IOException {
    BitSet b = squeezeMethod(remains, bytes);
    remains = writeSet(b);
  }

  public void finish() {}
}
