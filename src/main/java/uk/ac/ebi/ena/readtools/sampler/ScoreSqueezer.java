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

public class ScoreSqueezer {
  ScoreSqueezer() {}

  byte[] squeeze(byte value[]) {
    byte[] result = new byte[value.length / 2 + value.length % 2];
    int r = 0;
    for (int i = 0; i < value.length; i += 2)
      result[r++] = (byte) squeezePair(value[i], i + 1 < value.length ? value[i + 1] : (byte) 'A');

    return result;
  }

  int squeezePair(byte one, byte two) {
    return (translate(one) << 4) | (0xF & translate(two));
  }

  int translate(int value) {
    switch (value) {
      case 'A':
        return 0x00;
      case 'C':
        return 0x01;
      case 'G':
        return 0x02;
      case 'T':
        return 0x03;
      case 'N':
        return 0x04;

      default:
        throw new RuntimeException(String.format("byte %s could not be translated", value));
    }
  }
}
