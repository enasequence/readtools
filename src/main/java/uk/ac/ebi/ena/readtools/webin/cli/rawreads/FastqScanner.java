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

import java.io.File;
import java.io.InputStream;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.stream.Collectors;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import uk.ac.ebi.ena.readtools.loader.common.converter.AutoNormalizeQualityReadConverter;
import uk.ac.ebi.ena.readtools.loader.common.converter.ConverterException;
import uk.ac.ebi.ena.readtools.loader.fastq.FastqIterativeWriter;
import uk.ac.ebi.ena.readtools.loader.fastq.FastqIterativeWriter.READ_TYPE;
import uk.ac.ebi.ena.readtools.loader.fastq.PairedRead;
import uk.ac.ebi.ena.readtools.utils.Utils;
import uk.ac.ebi.ena.webin.cli.validator.message.ValidationMessage;
import uk.ac.ebi.ena.webin.cli.validator.message.ValidationOrigin;
import uk.ac.ebi.ena.webin.cli.validator.message.ValidationResult;


public abstract class
FastqScanner {
    private static final int MAX_LABEL_SET_SIZE = 10;

    private static final int PAIRING_THRESHOLD = 20;

    //TODO remove duplication
    private static final int DEFAULT_PRINT_FREQ = 1_000;
    private static final String PRINT_FREQ_PROPERTY_NAME = "webincli.scanner.print.freq";
    private static final int PRINT_FREQ = Integer.valueOf(System.getProperty(PRINT_FREQ_PROPERTY_NAME, String.valueOf(DEFAULT_PRINT_FREQ)));

    private static final Logger log = LoggerFactory.getLogger(FastqScanner.class);

    private final int expected_size;
    private final Set<String> labelSet = new HashSet<>();
    private final AtomicBoolean paired = new AtomicBoolean();

    private final Long readLimit;
    private Long count = 0L;

    abstract protected void logFlushMsg(String message);

    public FastqScanner(int expected_size) {
        this(null, expected_size);
    }

    /**
     * @param readLimit Only process limited number of reads per file. Prevents from processing all reads in a given file.
     */
    public FastqScanner(Long readLimit) {
        this(readLimit, readLimit.intValue() * 5);
    }

    /**
     * @param readLimit     Only process limited number of reads per file. Prevents from processing all reads in a given file.
     * @param expected_size
     */
    public FastqScanner(Long readLimit, int expected_size) {
        this.readLimit = readLimit;
        this.expected_size = expected_size;
    }

    public boolean
    getPaired() {
        return this.paired.get();
    }

    public void checkFiles(ValidationResult validationResult, RawReadsFile... rawReadsFiles) throws Throwable {
        if (null == rawReadsFiles || rawReadsFiles.length == 0) {
            //terminal error
            validationResult.add(ValidationMessage.error("No file provided."));
            return;
        }

        /**
         * Duplicate read name check across files has been dropped. Now, the check will be performed for read names
         * within the file only.
         */

        RawReadsFile mainFile = rawReadsFiles[0];

        /** Should ideally have high number of duplicates as it will point to higher pairing percentage.
         * To keep memory consumption lower, because we can tolerate false positive here, use lower expected read size. */
        BloomWrapper mainFileOnlyPairingBloomWrapper = new BloomWrapper(expected_size / 10);

        long mainFileReadCount =
                checkForDuplicatesAndAddPairingInformation(validationResult, mainFile, mainFileOnlyPairingBloomWrapper);

        if (!validationResult.isValid() || rawReadsFiles.length == 1) {
            return;
        }

        List<PairedFiles> pairedFiles = new ArrayList<>();

        //Note that i = 1. Meaning, we are now going to deal with every other file in the group as main file has already
        //been processed above.
        for (int i = 1; i < rawReadsFiles.length; i++) {
            RawReadsFile currentFile = rawReadsFiles[i];

            //Make a copy of the main file's pairing information so we do not have re-create it for every other file.
            BloomWrapper currentAndMainFilePairingBloomWrapper = mainFileOnlyPairingBloomWrapper.getCopy();

            //Check for duplicates in the current file and add its pairing information to main file's pairing information
            //to determine how well these two are paired.
            long currentFileReadCount = checkForDuplicatesAndAddPairingInformation(
                    validationResult, currentFile, currentAndMainFilePairingBloomWrapper);

            if (!validationResult.isValid()) {
                break;
            }

            //Calculate pairing percentage between current file and main file to determine how well they are paired.

            long readCount = Math.max(mainFileReadCount, currentFileReadCount);

            long pairedCount = currentAndMainFilePairingBloomWrapper.getPossibleDuplicateCount();

            double pairingPercentage = 100 * ((double) pairedCount / (double) readCount);

            pairedFiles.add(new PairedFiles(
                    mainFile.getFilename(), currentFile.getFilename(), pairingPercentage));
        }

        if (!validationResult.isValid())
            return;

        //Label set size and low pairing percentage validation.

        if (labelSet.size() <= rawReadsFiles.length) {
            paired.set(true);

            PairedFiles lowestPairingPercentagePair = pairedFiles.stream()
                    .sorted(Comparator.comparingDouble(pair -> pair.pairingPercentage))
                    .findFirst().get();

            //TODO: estimate bloom false positives impact
            if (lowestPairingPercentagePair.pairingPercentage < (double) PAIRING_THRESHOLD) {
                //terminal error
                validationResult.add(ValidationMessage.error(
                        String.format("Detected paired fastq submission with less than %d%% of paired reads between %s and %s",
                                PAIRING_THRESHOLD, lowestPairingPercentagePair.fileName1, lowestPairingPercentagePair.fileName2)));
            }
        } else if (labelSet.size() > rawReadsFiles.length) {
            validationResult.add(ValidationMessage.error(String.format(
                    "When submitting paired reads using two Fastq files the reads must follow Illumina paired read naming conventions. "
                            + "This was not the case for the submitted Fastq files: %s. Unable to determine pairing from set: %s",
                    rawReadsFiles,
                    labelSet.stream().limit(10).collect(Collectors.joining(",", "", 10 < labelSet.size() ? "..." : "")))));
        }
    }


    /**
     * Check the file for duplicate read names and updates the given bloom wrapper with pairing information.
     *
     * @param validationResult
     * @param rawReadsFile
     * @param pairingBloomWrapper
     * @return Number of read processed from the given file.
     * @throws Throwable
     */
    private long checkForDuplicatesAndAddPairingInformation(
            ValidationResult validationResult, RawReadsFile rawReadsFile, BloomWrapper pairingBloomWrapper) throws Throwable {

        ValidationResult fileValidationResult = rawReadsFile.getReportFile() == null
                ? validationResult.create(new ValidationOrigin("file", rawReadsFile.getFilename()))
                : validationResult.create(rawReadsFile.getReportFile().toFile(), new ValidationOrigin("file", rawReadsFile.getFilename()));

        /** Should ideally have a low to 0 number of duplicates. */
        BloomWrapper duplicationsBloomWrapper = new BloomWrapper(expected_size);

        Set<String> fileLabelSet = new HashSet<>();
        checkSingleFile(fileValidationResult, rawReadsFile, fileLabelSet, pairingBloomWrapper, duplicationsBloomWrapper);
        labelSet.addAll(fileLabelSet);

        //extra check for suspected reads
        if (fileValidationResult.isValid() && duplicationsBloomWrapper.hasPossibleDuplicates()) {
            // read name, list
            Map<String, Set<String>> duplicates = findAllduplications(duplicationsBloomWrapper, 100, rawReadsFile);

            ValidationResult duplicationResult = validationResult.create();
            duplicates.entrySet().stream().forEach(e -> duplicationResult.add(ValidationMessage.error(
                    String.format("Multiple (%d) occurrences of read name \"%s\" at: %s\n",
                            e.getValue().size(),
                            e.getKey(),
                            e.getValue().toString()))));

            if (duplicationResult.isValid()) {
                validationResult.add(ValidationMessage.info("No actual duplicate read names found."));
            }
        }

        return duplicationsBloomWrapper.getAddCount();
    }

    private void
    checkSingleFile(ValidationResult fileResult,
                    RawReadsFile rawReadsFile,
                    Set<String> labelSet,
                    BloomWrapper pairingBloomWrapper,
                    BloomWrapper duplicationsBloomWrapper) throws Throwable {
        try {
            read(rawReadsFile, labelSet, pairingBloomWrapper, duplicationsBloomWrapper);

            logFlushMsg("Processing completed. Result: OK\n");

            fileResult.add(ValidationMessage.info(String.format("Collected %d reads", count)));
            fileResult.add(ValidationMessage.info(String.format("Collected %d read labels: %s", labelSet.size(), labelSet)));
            fileResult.add(ValidationMessage.info(String.format("Has possible duplicate read name(s): " + duplicationsBloomWrapper.hasPossibleDuplicates())));
        } catch (ConverterException converterException) {
            converterException.printStackTrace();

            logFlushMsg(String.format("Processing completed. Result: %s\n", converterException));

            ValidationMessage dataProducerError = ValidationMessage.error(converterException.getMessage());
            dataProducerError.appendOrigin(new ValidationOrigin("line number", converterException.getLineNo()));
            fileResult.add(dataProducerError);
        }
    }

    private void
    read(RawReadsFile readsFile,
         Set<String> labels,
         BloomWrapper pairingBloomWrapper,
         BloomWrapper duplicationsBloomWrapper) throws Throwable {
        try (InputStream inputStream = Utils.openFastqInputStream(Paths.get(readsFile.getFilename()))) {
            String streamName = readsFile.getFilename();
            FastqReadScanner fastqReadScanner = new FastqReadScanner(
                    streamName, labels, pairingBloomWrapper, duplicationsBloomWrapper,
                    MAX_LABEL_SET_SIZE, PRINT_FREQ);
            AutoNormalizeQualityReadConverter readScanningConverter = new AutoNormalizeQualityReadConverter(
                    inputStream,
                    fastqReadScanner,
                    readLimit, "", readsFile.getFilename());

            log.info("Processing file " + readsFile.getFilename());
            readScanningConverter.run();

            count = readScanningConverter.getReadCount();
            logProcessedReadNumber(count);
            if (count <= 0) {
                throw new ConverterException( 0, "Empty file" );
            }
        }
    }

    protected abstract void logProcessedReadNumber(Long count);

    private Map<String, Set<String>> findAllduplications(BloomWrapper duplications, int limit, RawReadsFile rf) {
        Map<String, Integer> counts = new HashMap<>(limit);
        Map<String, Set<String>> results = new LinkedHashMap<>(limit);

        String msg = "Verifying possible duplicates for file " + rf.getFilename();
        log.info(msg);

        FastqIterativeWriter wrapper = new FastqIterativeWriter();
        wrapper.setFiles(new File[]{new File(rf.getFilename())});
        wrapper.setReadType(READ_TYPE.SINGLE);
        wrapper.setReadLimit(readLimit);

        Iterator<String> read_name_iterator = new DelegateIterator<PairedRead, String>(wrapper.iterator()) {
            @Override
            public String convert(PairedRead obj) {
                return obj.forward.getName();
            }
        };

        long index = 1;

        while (read_name_iterator.hasNext()) {
            String read_name = read_name_iterator.next();
            if (duplications.getPossibleDuplicates().contains(read_name)) {
                counts.put(read_name, counts.getOrDefault(read_name, 0) + 1);
                Set<String> dlist = results.getOrDefault(read_name, new LinkedHashSet<>());
                dlist.add(rf.getFilename() + ", read " + index);
                results.put(read_name, dlist);
            }
            ++index;
        }

        return results.entrySet()
                .stream()
                //only read names occurring more than once are considered duplicates
                .filter(e -> counts.get(e.getKey()) > 1)
                .limit(limit)
                .collect(Collectors.toMap(e -> e.getKey(), e -> e.getValue(), (v1, v2) -> v1, LinkedHashMap::new));
    }

    private static class PairedFiles {
        public String fileName1;
        public String fileName2;

        public double pairingPercentage;

        public PairedFiles(String fileName1, String fileName2, double pairingPercentage) {
            this.fileName1 = fileName1;
            this.fileName2 = fileName2;
            this.pairingPercentage = pairingPercentage;
        }
    }
}
