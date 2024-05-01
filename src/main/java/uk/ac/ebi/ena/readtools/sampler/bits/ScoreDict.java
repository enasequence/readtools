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

public class ScoreDict implements ConstantLengthDataDictionary<Character> {
  //    33[!]  0:     50037:   1.000000        0.500000        : **********
  //    34["]  1:         0:   0.794328        0.442688        :
  //    35[#]  2:   1480648:   0.630957        0.386863        : **************
  //    36[$]  3:       451:   0.501187        0.333861        : ******
  //    37[%]  4:      2138:   0.398107        0.284747        : *******
  //    38[&]  5:      4044:   0.316228        0.240253        : ********
  //    39[']  6:      6715:   0.251189        0.200760        : ********
  //    40[(]  7:     13140:   0.199526        0.166338        : *********
  //    41[)]  8:     17294:   0.158489        0.136807        : *********
  //    42[*]  9:     14716:   0.125893        0.111816        : *********
  //    43[+] 10:     16470:   0.100000        0.090909        : *********
  //    44[,] 11:     39925:   0.079433        0.073588        : **********
  //    45[-] 12:     75172:   0.063096        0.059351        : ***********
  //    46[.] 13:     12332:   0.050119        0.047727        : *********
  //    47[/] 14:     12180:   0.039811        0.038287        : *********
  //    48[0] 15:     20836:   0.031623        0.030653        : *********
  //    49[1] 16:     16051:   0.025119        0.024503        : *********
  //    50[2] 17:     26568:   0.019953        0.019562        : **********
  //    51[3] 18:     34084:   0.015849        0.015602        : **********
  //    52[4] 19:     44554:   0.012589        0.012433        : **********
  //    53[5] 20:    194873:   0.010000        0.009901        : ************
  //    54[6] 21:     45241:   0.007943        0.007881        : **********
  //    55[7] 22:     52031:   0.006310        0.006270        : **********
  //    56[8] 23:     49703:   0.005012        0.004987        : **********
  //    57[9] 24:     59774:   0.003981        0.003965        : **********
  //    58[:] 25:    346165:   0.003162        0.003152        : ************
  //    59[;] 26:     85035:   0.002512        0.002506        : ***********
  //    60[<] 27:    103149:   0.001995        0.001991        : ***********
  //    61[=] 28:    467880:   0.001585        0.001582        : *************
  //    62[>] 29:    173698:   0.001259        0.001257        : ************
  //    63[?] 30:    895855:   0.001000        0.000999        : *************
  //    64[@] 31:    541613:   0.000794        0.000794        : *************
  //    65[A] 32:    962739:   0.000631        0.000631        : *************
  //    66[B] 33:   1454517:   0.000501        0.000501        : **************
  //    67[C] 34:   1251956:   0.000398        0.000398        : **************
  //    68[D] 35:   1900604:   0.000316        0.000316        : **************
  //    69[E] 36:   5344916:   0.000251        0.000251        : ***************
  //    70[F] 37:   3118306:   0.000200        0.000199        : **************
  //    71[G] 38:  11949190:   0.000158        0.000158        : ****************

  @Override
  public int getIncomingBlockLengthBytes() {
    return 1;
  }

  @Override
  public int getTranslatedLength() {
    return 6;
  }

  @Override
  public int translate(int value) {
    if (33 <= value && 126 >= value) {
      switch (value) {
        case 71:
          return 0; // 71[G]:  11949190
        case 69:
          return 1; // 69[E]:   5344916
        case 70:
          return 2; // 70[F]:   3118306
        case 68:
          return 3; // 68[D]:   1900604
        case 35:
          return 4; // 35[#]:   1480648
        case 66:
          return 5; // 66[B]:   1454517
        case 67:
          return 6; // 67[C]:   1251956
        case 65:
          return 7; // 65[A]:    962739
        case 63:
          return 8; // 63[?]:    895855
        case 64:
          return 9; // 64[@]:    541613
        case 61:
          return 10; // 61[=]:    467880
        case 58:
          return 11; // 58[:]:    346165
        case 53:
          return 12; // 53[5]:    194873
        case 62:
          return 13; // 62[>]:    173698
        case 60:
          return 14; // 60[<]:    103149
        case 59:
          return 15; // 59[;]:     85035
        case 45:
          return 16; // 45[-]:     75172
        case 57:
          return 17; // 57[9]:     59774
        case 55:
          return 18; // 55[7]:     52031
        case 33:
          return 19; // 33[!]:     50037
        case 56:
          return 20; // 56[8]:     49703
        case 54:
          return 21; // 54[6]:     45241
        case 52:
          return 22; // 52[4]:     44554
        case 44:
          return 23; // 44[,]:     39925
        case 51:
          return 24; // 51[3]:     34084
        case 50:
          return 25; // 50[2]:     26568
        case 48:
          return 26; // 48[0]:     20836
        case 41:
          return 27; // 41[)]:     17294
        case 43:
          return 28; // 43[+]:     16470
        case 49:
          return 29; // 49[1]:     16051
        case 42:
          return 30; // 42[*]:     14716
        case 40:
          return 31; // 40[(]:     13140
        case 46:
          return 32; // 46[.]:     12332
        case 47:
          return 33; // 47[/]:     12180
        case 39:
          return 34; // 39[']:      6715
        case 38:
          return 35; // 38[&]:      4044
        case 37:
          return 36; // 37[%]:      2138
        case 36:
          return 37; // 36[$]:       451
        case 34:
          return 38; // 34["]:         0
        default:
          return 38;
      }
    } else throw new RuntimeException();
  }

  @Override
  public int translateBack(int value) {
    if (33 <= value && 126 >= value) {
      return value;

    } else throw new RuntimeException();
  }

  @Override
  public int getDictionarySize() {

    return 5;
  }
}
