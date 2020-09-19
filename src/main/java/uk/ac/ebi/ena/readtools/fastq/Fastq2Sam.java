/*
 * Copyright 2010-2020 EMBL - European Bioinformatics Institute
 * Licensed under the Apache License, Version 2.0 (the "License"); you may not use this
 * file except in compliance with the License. You may obtain a copy of the License at
 * http://www.apache.org/licenses/LICENSE-2.0
 * Unless required by applicable law or agreed to in writing, software distributed under the
 * License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
 * CONDITIONS OF ANY KIND, either express or implied. See the License for the
 * specific language governing permissions and limitations under the License.
 */
package uk.ac.ebi.ena.readtools.fastq;

/**
 * A wrapper class around picard.sam.FastqToSam.
 */
public class Fastq2Sam {

    private final String sampleName;

    private final String inputFile1;
    private final String inputFile2;

    private final String outputFile;

    private final String tempDir;

    /**
     *
     * @param sampleName - Must not be null or empty. Used as the value of SAM header 'SM'.
     * @param inputFile - Path to input Fastq file. File compression is determined automatically.
     * @param outputFile - Path to output SAM file. File extension determines SAM format.
     * @param tempDir - Path to the temporary directory to be used during the process. Pass null to use system default.
     */
    public Fastq2Sam(String sampleName, String inputFile, String outputFile, String tempDir) {
        this(sampleName, inputFile, null, outputFile, tempDir);
    }

    /**
     * For paired Fastq files.
     *
     * @param sampleName - Must not be null or empty. Used as the value of SAM header 'SM'
     * @param inputFile1 - Path to input Fastq file. File compression is determined automatically.
     * @param inputFile2 - Path to input Fastq file. File compression is determined automatically.
     * @param outputFile - Path to output SAM file.  File extension determines SAM format.
     * @param tempDir - Path to the temporary directory to be used during the process.. Pass null to use system default.
     */
    public Fastq2Sam(String sampleName, String inputFile1, String inputFile2, String outputFile, String tempDir) {
        this.sampleName = sampleName;
        this.inputFile1 = inputFile1;
        this.inputFile2 = inputFile2;
        this.outputFile = outputFile;
        this.tempDir = tempDir;

        if (sampleName == null || sampleName.trim().isEmpty()) {
            throw new IllegalArgumentException("Sample name must not be null or empty.");
        }
    }

    public void convert() {
        new picard.sam.FastqToSam().instanceMain(new String[]{
                "SM=" + sampleName,
                "F1=" + inputFile1,
                "F2=" + inputFile2,
                "O=" + outputFile,
                "TMP_DIR=" + tempDir,
                "ALLOW_AND_IGNORE_EMPTY_LINES=true"
        });
    }
}
