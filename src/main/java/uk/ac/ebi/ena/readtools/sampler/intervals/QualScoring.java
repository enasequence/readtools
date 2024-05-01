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
package uk.ac.ebi.ena.readtools.sampler.intervals;

import java.util.HashMap;
import java.util.Map;

public enum QualScoring {
  PHRED(1),
  LOGODDS(2);

  static final double delta = 1e-8;
  static final int hbound = 41;
  static final int lbound = -5;
  final double[] probs = new double[hbound - lbound];
  double math_e;
  double disp;

  private int getPhredQ(double p) {
    return (int) (-10 * Math.log10(p));
  }

  private int getLogOddsQ(double p) {
    return (int) (-10 * Math.log10(p / (1 - p)));
  }

  private QualScoring(int type) {

    int count = 0;
    int ps[] = new int[hbound - lbound];

    switch (type) {
      case 1:
        for (double p = 0; p < 1; p = p + delta) {
          int value = getPhredQ(p);
          if (lbound <= value && hbound > value) {
            ps[value - lbound] += 1;
            count++;
          }
        }
        break;
      case 2:
        for (double p = 0; p < 1; p = p + delta) {
          int value = getLogOddsQ(p);
          if (lbound <= value && hbound > value) {
            ps[value - lbound] += 1;
            count++;
          }
        }
        break;
      default:
        throw new RuntimeException();
    }

    double sum = 0;
    for (int i = 0; i < ps.length; ++i) {
      probs[i] = count == 0 ? 0 : (double) ps[i] / (double) count;
      sum += probs[i];
      math_e += probs[i] * (i + lbound);
      //      System.out.printf( "%d:probs[ %d ] = %f\n", type, i + lbound, probs[ i ] );
    }

    for (int i = 0; i < ps.length; ++i) disp += Math.pow((i + lbound) - math_e, 2) * probs[i];
    disp = Math.sqrt(disp);
    System.out.printf(
        "sum: %f, exp: %f, disp: %f, avg: %f\n", sum, math_e, disp, sum / (double) ps.length);
  }

  public double calcWeight(long[] row, int offset) {
    double result = 0;
    long tcount = 0;

    for (int i = 0; i < probs.length; ++i) {
      int index = offset + i + lbound;
      if (index < 0) // || ( ( i + lbound ) == 2  ) )  // 2 - Illumina 1.5 /B/
      continue;
      tcount += row[index];
    }

    /*
            for( int i = 0; i < probs.length; ++i )
            {
                int index = offset + i + lbound;
                if( index < 0 )
                    continue;
                row_norm = (double)row[ index ] / (double)tcount;
            }
    */
    int scount = 0;
    long count = 0;
    double e = 0;
    double d = 0;
    for (int i = 0; i < probs.length; ++i) {
      int index = offset + i + lbound;
      if (index < 0) count = 0;
      else count = row[index];

      e += (i + lbound) * ((double) count / (double) tcount);
    }

    for (int i = 0; i < probs.length; ++i) {
      int index = offset + i + lbound;

      if (index < 0)
        //                || ( i + lbound ) == 2 )// 2 - Illumina 1.5 /B/
        count = 0;
      else count = row[index];

      d += Math.pow(i + lbound - math_e, 2) * ((double) count / (double) tcount);

      double p = (double) count / (double) tcount;
      double sd = /*count == 0 ? 0 :*/ Math.pow((probs[i] - p), 2);

      System.out.printf(
          "%7s: probs[%3d] = %f\t%f\t%f\n", this.toString(), i + lbound, probs[i], p, sd);

      result += sd;
      scount++;
    }

    System.out.printf("exp: %f\tdisp: %f\n", e, Math.sqrt(d));

    // return Math.sqrt( d );
    return Math.sqrt((double) result / (double) scount);
  }

  public static Map<QualScoring, Double> forRow(long[] row, int offset) {
    Map<QualScoring, Double> map = new HashMap<QualScoring, Double>();
    for (QualScoring qs : QualScoring.values()) map.put(qs, qs.calcWeight(row, offset));

    return map;
  }
}
