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
package uk.ac.ebi.ena.readtools.webin.cli.rawreads;

import uk.ac.ebi.ena.readtools.loader.common.writer.ReadWriter;
import uk.ac.ebi.ena.readtools.loader.common.writer.ReadWriterException;
import uk.ac.ebi.ena.readtools.loader.common.writer.Spot;
import uk.ac.ebi.ena.readtools.loader.fastq.Read;

public abstract class InsdcStandardCheckingScanner implements ReadWriter<Read, Spot> {
  /*
  A	Adenine
  C	Cytosine
  G	Guanine
  T (or U)	Thymine (or Uracil)
  R	A or G
  Y	C or T
  S	G or C
  W	A or T
  K	G or T
  M	A or C
  B	C or G or T
  D	A or G or T
  H	A or C or T
  V	A or C or G
  N	any base
  . or -
  */
  private static final String IUPAC_CODES = "ACGTURYSWKMBDHVN.-";

  private final int printFreq;
  private Long count = 0L;

  protected abstract void logProcessedReadNumber(long cnt);

  public InsdcStandardCheckingScanner(int printFreq) {
    this.printFreq = printFreq;
  }

  @Override
  public void cascadeErrors() throws ReadWriterException {
    // Not implemented
  }

  @Override
  public void write(Read read) throws ReadWriterException {
    // Check for empty reads
    //        if (read == null || read.getBases() == null || read.getBases().isEmpty() ||
    // read.getBases().trim().isEmpty()) {
    //            throw new ReadWriterException("File must not contain any empty reads.");
    //        }
    //
    //        // Check for valid IUPAC codes and no more than 50% non-AUTCG bases
    //        long validBasesCount = read.getBases().chars().filter(c -> IUPAC_CODES.indexOf(c) !=
    // -1).count();
    //        long nonAUTCGCount = read.getBases().chars().filter(c -> "RYKMSWBDHVN".indexOf(c) !=
    // -1).count();
    //
    //        if (validBasesCount < read.getBaseCount() || nonAUTCGCount > read.getBaseCount() / 2)
    // {
    //            throw new ReadWriterException("Reads must contain only valid IUPAC codes with no
    // more than 50% non-AUTCG bases.");
    //        }
    //
    //        // When file contains base quality scores
    //        if (read.getQualityScores() != null && !read.getQualityScores().isEmpty()) {
    //            // Compute the average quality score for the read
    //            double averageQuality = read.getQualityScores().chars()
    //                    .mapToDouble(c -> c - '!').average().orElse(0.0);
    //
    //            // Check if the average quality is below the threshold
    //            if (averageQuality < 30) {
    //                throw new ReadWriterException("At least 50% of reads must have an average
    // quality score >= 30.");
    //            }
    //        }

    if (read.getBases().length() != read.getQualityScores().length()) {
      throw new ReadWriterException(
          "Mismatch between length of read bases and qualities",
          ReadWriterException.ErrorType.BASES_QUALITIES_LENGTH_MISMATCH);
    }

    count++;
    if (0 == count % printFreq) {
      logProcessedReadNumber(count);
    }
  }

  @Override
  public void setWriter(ReadWriter readWriter) {
    throw new RuntimeException("Not implemented");
  }
}
