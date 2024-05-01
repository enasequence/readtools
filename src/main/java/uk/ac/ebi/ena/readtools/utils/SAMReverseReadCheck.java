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
package uk.ac.ebi.ena.readtools.utils;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import java.io.File;
import uk.ac.ebi.ena.readtools.sam.Sam2Fastq;

public class SAMReverseReadCheck {

  private final File samFile;

  public static void main(String[] args) {
    if (args.length == 0 || args[0].isEmpty()) {
      System.out.println("Provide path to a SAM/BAM/CRAM file.");
      System.exit(1);
    }

    File samFile = new File(args[0]);
    if (!samFile.exists() || samFile.length() == 0) {
      System.out.println("Input SAM/BAM/CRAM file does not exists or is empty.");
      System.exit(1);
    }

    Output output = new SAMReverseReadCheck(samFile).check();
    if (output == null) {
      System.out.println("No reverse reads found.");
    } else {
      System.out.println(
          String.format(
              "%s\n%s\n%s\n%d", output.readName, output.bases, output.qualities, output.flags));
    }
  }

  public SAMReverseReadCheck(File samFile) {
    this.samFile = samFile;
  }

  public Output check() {
    Output result = null;

    SamReaderFactory.setDefaultValidationStringency(ValidationStringency.LENIENT);

    try (SamReader samReader = SamReaderFactory.makeDefault().open(samFile)) {
      for (final SAMRecord currentRecord : samReader) {

        if (currentRecord.isSecondaryOrSupplementary()
            && !Sam2Fastq.INCLUDE_NON_PRIMARY_ALIGNMENTS) {
          continue;
        }

        // Skip non-PF reads as necessary
        if (currentRecord.getReadFailsVendorQualityCheckFlag() && !Sam2Fastq.INCLUDE_NON_PF_READS) {
          continue;
        }

        if (currentRecord.getReadNegativeStrandFlag()) {
          result = new Output();
          result.flags = currentRecord.getFlags();
          result.readName = currentRecord.getReadName();
          result.bases = currentRecord.getReadString();
          result.qualities = currentRecord.getBaseQualityString();
          break;
        }
      }
    } catch (Exception ex) {
      throw new RuntimeException(ex);
    }

    return result;
  }

  class Output {

    public int flags;
    public String readName, bases, qualities;
  }
}
