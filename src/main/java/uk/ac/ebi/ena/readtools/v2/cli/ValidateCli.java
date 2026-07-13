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
package uk.ac.ebi.ena.readtools.v2.cli;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import java.io.File;
import java.util.List;
import java.util.stream.Collectors;
import uk.ac.ebi.ena.readtools.v2.FileFormat;
import uk.ac.ebi.ena.readtools.v2.validator.ReadsValidationException;
import uk.ac.ebi.ena.readtools.v2.validator.ValidatorWrapper;

/**
 * Standalone CLI wrapper around {@link ValidatorWrapper}, running the same client-side read
 * validation that webin-cli performs, without the manifest / upload machinery.
 *
 * <p>Usage: java -cp readtools-all.jar uk.ac.ebi.ena.readtools.v2.cli.ValidateCli \ --format FASTQ
 * [--full] file1 [file2]
 */
public class ValidateCli {
  // webin-cli's own limits: quick = first 100k reads, extended = first 100M reads.
  private static final long QUICK_READ_LIMIT = 100_000L;
  private static final long EXTENDED_READ_LIMIT = 100_000_000L;

  @Parameter(
      names = {"--format", "-f"},
      description = "File format: FASTQ, BAM or CRAM",
      required = true)
  private FileFormat format;

  @Parameter(
      names = {"--full"},
      description =
          "Validate up to 100M reads (webin-cli 'extended' mode). Default is quick: first 100k.")
  private boolean full = false;

  @Parameter(description = "<file> [file2]  (two files = paired FASTQ)", required = true)
  private List<String> files = new java.util.ArrayList<>();

  @Parameter(
      names = {"--help", "-h"},
      help = true)
  private boolean help = false;

  public static void main(String[] args) {
    ValidateCli cli = new ValidateCli();
    JCommander jc = JCommander.newBuilder().addObject(cli).programName("readtools-validate").build();

    try {
      jc.parse(args);
    } catch (Exception e) {
      System.err.println("ERROR: " + e.getMessage());
      jc.usage();
      System.exit(2);
    }

    if (cli.help) {
      jc.usage();
      return;
    }

    System.exit(cli.run());
  }

  private int run() {
    List<File> fileList = files.stream().map(File::new).collect(Collectors.toList());
    for (File f : fileList) {
      if (!f.isFile()) {
        System.err.println("ERROR: file not found: " + f.getPath());
        return 2;
      }
    }

    long limit = full ? EXTENDED_READ_LIMIT : QUICK_READ_LIMIT;
    ValidatorWrapper wrapper = new ValidatorWrapper(fileList, format, limit);

    try {
      wrapper.run();
    } catch (ReadsValidationException e) {
      System.out.println("RESULT: INVALID");
      System.out.println("  " + e.getMessage());
      return 1;
    } catch (RuntimeException e) {
      System.out.println("RESULT: INVALID (file structure / parse error)");
      System.out.println("  " + (e.getCause() != null ? e.getCause().getMessage() : e.getMessage()));
      return 1;
    }

    System.out.println("RESULT: VALID");
    System.out.println("  format:       " + format);
    System.out.println("  mode:         " + (full ? "full (<=100M reads)" : "quick (<=100k reads)"));
    if (format == FileFormat.FASTQ && fileList.size() > 1) {
      System.out.println("  paired:       " + wrapper.isPaired());
    }
    wrapper
        .getFileQualityStats()
        .forEach(
            s ->
                System.out.printf(
                    "  %s: reads=%d, highQuality(avg>=30)=%d%n",
                    s.getFile().getName(), s.getReadCount(), s.getHighQualityReadCount()));
    return 0;
  }
}
