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

public class BaseDict implements ConstantLengthDataDictionary<Character> {

  @Override
  public int getIncomingBlockLengthBytes() {
    return 1;
  }

  @Override
  public int getTranslatedLength() {
    return 3;
  }

  /*
      7   111
      a   4   010
      e   4   000
      f   3   1101
      h   2   1010

      i   2   1000
      m   2   0111
      n   2   0010
      s   2   1011
      t   2   0110
  */
  @Override
  public int translate(int value) {
    switch (value) {
      case 'A':
        return 0x1;
      case 'T':
        return 0x3;
      case 'G':
        return 0x5;
      case 'C':
        return 0x7;
      case 'N':
        return 0x0;

      default:
        throw new RuntimeException(String.format("Unknown value: 0x%h", value));
    }
  }

  @Override
  public int translateBack(int value) {
    switch (value) {
      case 0x1:
        return 'A';
      case 0x3:
        return 'T';
      case 0x5:
        return 'G';
      case 0x7:
        return 'C';
      case 0x0:
        return 'N';
      default:
        throw new RuntimeException(String.format("Unknown value: 0x%h", value));
    }
  }

  @Override
  public int getDictionarySize() {

    return 5;
  }
}
