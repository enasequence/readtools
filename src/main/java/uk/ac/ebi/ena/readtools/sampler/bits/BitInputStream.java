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
package uk.ac.ebi.ena.readtools.sampler.bits;

import java.io.EOFException;
import java.io.FilterInputStream;
import java.io.IOException;
import java.io.InputStream;

public class BitInputStream extends FilterInputStream {
  final int BYTE_LENGTH = 8;
  int rlength = 0;
  int rvalue = 0;

  static final int[] bit_mask =
      new int[] {
        ~(0xFFFFFFFF << 0), ~(0xFFFFFFFF << 1), ~(0xFFFFFFFF << 2), ~(0xFFFFFFFF << 3),
        ~(0xFFFFFFFF << 4), ~(0xFFFFFFFF << 5), ~(0xFFFFFFFF << 6), ~(0xFFFFFFFF << 7),
        ~(0xFFFFFFFF << 8), ~(0xFFFFFFFF << 9), ~(0xFFFFFFFF << 10), ~(0xFFFFFFFF << 11),
        ~(0xFFFFFFFF << 12), ~(0xFFFFFFFF << 13), ~(0xFFFFFFFF << 14), ~(0xFFFFFFFF << 15),
        ~(0xFFFFFFFF << 16), ~(0xFFFFFFFF << 17), ~(0xFFFFFFFF << 18), ~(0xFFFFFFFF << 19),
        ~(0xFFFFFFFF << 20), ~(0xFFFFFFFF << 21), ~(0xFFFFFFFF << 22), ~(0xFFFFFFFF << 23),
        ~(0xFFFFFFFF << 24), ~(0xFFFFFFFF << 25), ~(0xFFFFFFFF << 26), ~(0xFFFFFFFF << 27),
        ~(0xFFFFFFFF << 28), ~(0xFFFFFFFF << 29), ~(0xFFFFFFFF << 30), ~(0xFFFFFFFF << 31)
      };

  public BitInputStream(InputStream stream) {
    super(stream);
  }

  private int readbyte() throws IOException {
    int r = super.read();
    if (-1 == r) throw new EOFException();
    return r & bit_mask[BYTE_LENGTH];
  }

  public int read(int length) throws IOException {
    int upto = (length - rlength) / BYTE_LENGTH;

    int bits = rvalue;
    int blength = rlength;

    for (int i = 0; i < upto; ++i) {
      bits = bits | (readbyte() << blength);
      blength += BYTE_LENGTH;
    }

    rlength = blength - length;

    if (rlength >= 0) {
      rvalue = bits >> length;
    } else {
      int b = readbyte();
      rlength += BYTE_LENGTH;
      rvalue = b >> (BYTE_LENGTH - rlength);
      bits = bits | ((b & bit_mask[BYTE_LENGTH - rlength]) << (blength));
    }
    return bits & bit_mask[length];
  }
}
