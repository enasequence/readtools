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

import htsjdk.samtools.SAMException;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import uk.ac.ebi.ena.readtools.common.reads.CasavaRead;
import uk.ac.ebi.ena.readtools.v2.provider.FastqReadsProvider;
import uk.ac.ebi.ena.readtools.v2.provider.ReadsProviderFactory;
import uk.ac.ebi.ena.readtools.v2.read.FastqRead;
import uk.ac.ebi.ena.readtools.webin.cli.rawreads.BloomWrapper;

public class FastqReadsValidator extends ReadsValidator {
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
  SPACE HERE
  8    <read> Numerical Read number. 1 can be single read or read 2 of pairedend
  9    <is
      filtered>
      Y or N Y if the read is filtered, N otherwise
  10    <control
      number>
      Numerical 0 when none of the control bits are on, otherwise it is
      an even number. See below.
  11    <index
      sequence>
      ACTG Index sequence
      */
  //                                                          1        :  2   :    3       :   4  :
  // 5   :   6   :  7          8 :  9 :  10         : 11
  //    final static private Pattern p_casava_1_8_name = Pattern.compile(
  // "^@([a-zA-Z0-9_-]+:[0-9]+:[a-zA-Z0-9]+:[0-9]+:[0-9]+:[0-9-]+:[0-9-]+)
  // ([12]):[YN]:[0-9]*[02468]:[ACGTN]+$" );
  // relaxed regular expression — canonical definition in CasavaRead
  public static final Pattern P_CASAVA_18_NAME = CasavaRead.P_CASAVA_18_NAME;
  private static final Pattern pQuals = Pattern.compile("^([!-~]*?)$"); // qualities
  private ReadStyle readStyle = null; // Field to keep track of the read style

  public FastqReadsValidator(long readCountLimit) {
    super(readCountLimit);
  }

  @Override
  public boolean validate(ReadsProviderFactory readsProviderFactory)
      throws ReadsValidationException {
    BloomWrapper duplicationsBloomWrapper = new BloomWrapper(5 * readCountLimit);
    Map<String, Integer> counts = new HashMap<>(100);

    long readCount = 0;
    try {
      try (FastqReadsProvider provider =
          (FastqReadsProvider) readsProviderFactory.makeReadsProvider()) {
        for (FastqRead read : provider) {
          if (readCount >= readCountLimit) {
            break;
          }
          readCount++;

          if (readCount == 1) {
            determineReadStyle(read.getName()); // Determine style based on the first read
          }

          validateRead(read, readCount);

          duplicationsBloomWrapper.add(read.getName());

          extraReadsValidation(readStyle, readCount, read);
        }
      }

      try (FastqReadsProvider provider =
          (FastqReadsProvider) readsProviderFactory.makeReadsProvider()) {
        long dupCheckReadCount = 0;
        for (FastqRead read : provider) {
          if (dupCheckReadCount >= readCountLimit) {
            break;
          }
          dupCheckReadCount++;

          duplicationsBloomWrapper.add(read.getName());
          if (duplicationsBloomWrapper.getPossibleDuplicates().contains(read.getName())) {
            counts.put(read.getName(), counts.getOrDefault(read.getName(), 0) + 1);
          }
        }
      }

      if (duplicationsBloomWrapper.hasPossibleDuplicates()) {
        StringBuilder errorReport = new StringBuilder();
        boolean found = false;

        for (Map.Entry<String, Integer> e : counts.entrySet()) {
          if (e.getValue() > 1) {
            found = true;
            errorReport.append(
                String.format(
                    "Multiple (%d) occurrences of read name \"%s\"", e.getValue(), e.getKey()));
          }
        }

        if (found) {
          throw new ReadsValidationException(errorReport.toString());
        }
      }

      return true;
    } catch (SAMException e) {
      throw new ReadsValidationException(e.getMessage(), readCount);
    } catch (ReadsValidationException rve) {
      throw rve;
    } catch (Exception e) {
      throw new RuntimeException(e);
    }
  }

  protected void extraReadsValidation(ReadStyle readStyle, long readCount, FastqRead read)
      throws ReadsValidationException {}

  private void determineReadStyle(String name) {
    Matcher casavaMatcher = P_CASAVA_18_NAME.matcher(name);
    readStyle = casavaMatcher.matches() ? ReadStyle.CASAVA18 : ReadStyle.FASTQ;
  }

  private void validateRead(FastqRead read, long readCount) throws ReadsValidationException {
    validateReadName(read.getName(), readCount);
    validateQualityScores(read.getName(), read.getQualityScores(), readCount);
  }

  private void validateReadName(String name, long readCount) throws ReadsValidationException {
    switch (readStyle) {
      case CASAVA18:
        Matcher casavaMatcher = P_CASAVA_18_NAME.matcher(name);
        if (!casavaMatcher.matches()) {
          throw new ReadsValidationException("Invalid CASAVA 1.8 read name", readCount, name);
        }
        break;
      case FASTQ:
        // For FASTQ, no specific pattern, but you can add basic checks if needed
        if (name == null || name.trim().isEmpty()) {
          throw new ReadsValidationException("Invalid FASTQ read name", readCount, name);
        }
        // Other basic validations for FASTQ can be added here if needed
        break;
      default:
        throw new ReadsValidationException(
            "Read style not determined for read name", readCount, name);
    }
  }

  private void validateQualityScores(String readName, String qualityScores, long readCount)
      throws ReadsValidationException {
    Matcher matcher = pQuals.matcher(qualityScores);
    if (!matcher.matches()) {
      throw new ReadsValidationException(
          "Invalid quality scores: " + qualityScores, readCount, readName);
    }
  }

  public enum ReadStyle {
    FASTQ,
    CASAVA18
  }
}
