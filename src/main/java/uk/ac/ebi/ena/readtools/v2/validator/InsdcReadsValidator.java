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
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;
import org.apache.commons.lang3.StringUtils;
import uk.ac.ebi.ena.readtools.v2.FileFormat;
import uk.ac.ebi.ena.readtools.v2.provider.ReadsProvider;
import uk.ac.ebi.ena.readtools.v2.provider.ReadsProviderFactory;
import uk.ac.ebi.ena.readtools.v2.read.IRead;

public class InsdcReadsValidator extends ReadsValidator {
  public static final String IUPAC_CODES = "ACGTURYSWKMBDHVN.-";
  private final Set<Character> iupacSet;
  private static final int MIN_QUALITY_SCORE = 30;

  public static String ERROR_NULL_READS = "Reads cannot be null";
  public static String ERROR_NO_READS = "Submitted files must contain a minimum of 1 sequence read";
  public static String ERROR_EMPTY_READ = "Submitted files must not contain any empty reads";
  public static String ERROR_READ_NAME_LENGTH = "Read name length exceeds 256 characters";
  public static String ERROR_NOT_IUPAC = "Reads must contain only valid IUPAC codes";
  public static String ERROR_NOT_AUTCG =
      "Reads must contain only valid IUPAC codes, "
          + "with no more than 50% of bases being non-AUTCG";
  public static String ERROR_EMPTY_QUALITY =
      "Submitted files must not contain any empty quality score lines";
  public static String ERROR_BASES_QUALITIES_LENGTH_MISMATCH =
      "Mismatch between length of read bases and qualities";
  public static String ERROR_QUALITY =
      "When submitted file contains base quality scores "
          + "then >= 50% of reads must have average quality >= 30";
  public static String INVALID_FILE = "Invalid file structure";

  private long readCount = 0;
  private long highQualityReadCount = 0;

  public InsdcReadsValidator(long readCountLimit) {
    super(readCountLimit);

    iupacSet = new HashSet<>();
    for (char c : IUPAC_CODES.toCharArray()) {
      iupacSet.add(c);
    }
  }

  public long getReadCount() {
    return readCount;
  }

  public long getHighQualityReadCount() {
    return highQualityReadCount;
  }

  @Override
  public boolean validate(ReadsProviderFactory readsProviderFactory)
      throws ReadsValidationException {
    long autcgCount = 0;
    long basesCount = 0;
    FileFormat inputFormat = readsProviderFactory.getFormat();
    boolean samLikeFormat = inputFormat == FileFormat.BAM || inputFormat == FileFormat.CRAM;

    try (ReadsProvider<? extends IRead> provider = readsProviderFactory.makeReadsProvider()) {

      if (provider == null) {
        throw new ReadsValidationException(ERROR_NULL_READS);
      }

      Iterator<? extends IRead> iterator;
      try {
        iterator = provider.iterator();
      } catch (SAMException e) {
        throw new ReadsValidationException(INVALID_FILE + ": " + e.getMessage());
      }

      if (!iterator.hasNext()) {
        throw new ReadsValidationException(ERROR_NO_READS);
      }

      while (iterator.hasNext()) {
        if (readCount >= readCountLimit) {
          break;
        }

        IRead read = iterator.next();
        String bases = read.getBases();
        String qualityScores = read.getQualityScores();

        readCount++;

        String effectiveBases = bases;
        // In SAM/BAM/CRAM, "*" means sequence is absent (length 0), not a literal base character.
        if (samLikeFormat && "*".equals(bases)) {
          effectiveBases = "";
        }

        if (StringUtils.isBlank(effectiveBases)) {
          throw new ReadsValidationException(ERROR_EMPTY_READ, readCount);
        }
        if (qualityScores == null || StringUtils.isBlank(qualityScores)) {
          throw new ReadsValidationException(ERROR_EMPTY_QUALITY, readCount, read.getName());
        }

        String effectiveQualityScores = qualityScores;
        // In SAM/BAM/CRAM, "*" means quality is absent (length 0), not a literal quality character.
        if (samLikeFormat && "*".equals(qualityScores)) {
          effectiveQualityScores = "";
        }

        if (effectiveBases.length() != effectiveQualityScores.length()) {
          throw new ReadsValidationException(
              ERROR_BASES_QUALITIES_LENGTH_MISMATCH, readCount, read.getName());
        }

        if (read.getName().trim().length() > 256) {
          throw new ReadsValidationException(ERROR_READ_NAME_LENGTH, readCount, read.getName());
        }

        basesCount += effectiveBases.length();
        for (char base : effectiveBases.toUpperCase().toCharArray()) {
          if (iupacSet.contains(base)) {
            if (base == 'A' || base == 'U' || base == 'T' || base == 'C' || base == 'G') {
              autcgCount++;
            }
          } else {
            throw new ReadsValidationException(ERROR_NOT_IUPAC, readCount, effectiveBases);
          }
        }

        int totalQuality = 0;
        for (char q : effectiveQualityScores.toCharArray()) {
          totalQuality += q - '!'; // Phred+33 0 at !
        }
        if ((double) totalQuality / effectiveQualityScores.length() >= MIN_QUALITY_SCORE) {
          highQualityReadCount++;
        }
      }

      if ((basesCount - autcgCount) > (basesCount / 2)) {
        throw new ReadsValidationException(ERROR_NOT_AUTCG, readCount);
      }

      //      if ((double) highQualityReadCount / readCount < 0.5) {
      //        throw new ReadsValidationException(ERROR_QUALITY, readCount);
      //      }
    } catch (ReadsValidationException rve) {
      throw rve;
    } catch (Exception e) {
      throw new RuntimeException(e);
    }

    return true;
  }
}
