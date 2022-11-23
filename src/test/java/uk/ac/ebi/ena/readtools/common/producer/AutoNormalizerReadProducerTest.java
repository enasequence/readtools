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
package uk.ac.ebi.ena.readtools.common.producer;

import java.io.InputStream;
import java.nio.file.Path;
import java.nio.file.Paths;

import org.junit.Test;

import uk.ac.ebi.ena.readtools.loader.common.converter.AutoNormalizeQualityReadConverter;
import uk.ac.ebi.ena.readtools.loader.common.converter.ConverterException;
import uk.ac.ebi.ena.readtools.loader.common.writer.ReadWriter;
import uk.ac.ebi.ena.readtools.loader.common.writer.ReadWriterException;
import uk.ac.ebi.ena.readtools.loader.common.writer.Spot;
import uk.ac.ebi.ena.readtools.loader.fastq.Read;
import uk.ac.ebi.ena.readtools.utils.Utils;

public class AutoNormalizerReadProducerTest {
    private static final long FASTQ_VALIDATION_MAX_DURATION_MS = 4_000;

    @Test (timeout = FASTQ_VALIDATION_MAX_DURATION_MS + 1_000)
    public void test() throws Exception {
        Path filePath = Paths.get("src/test/resources/tst.fastq.bz2");

        try (InputStream is = Utils.openFastqInputStream(filePath)) {
            AutoNormalizeQualityReadConverter dp =
                    new AutoNormalizeQualityReadConverter(is, "", filePath.toString());
            dp.setWriter(new ReadWriter<Read, Spot>() {
                @Override
                public void cascadeErrors() throws ReadWriterException {
                }

                @Override
                public void
                write(Read spot) {
                }

                @Override
                public void setWriter(ReadWriter<Spot, ?> readWriter) {
                    throw new RuntimeException("Not implemented");
                }

                @Override
                public boolean isOk() {
                    return true;
                }
            });
            dp.run();
            if (!dp.isOk()) {
                if (dp.getStoredException() instanceof ConverterException) {
                    dp.getStoredException().printStackTrace();
                    throw new Exception("DataProducerException");
                } else {
                    dp.getStoredException().printStackTrace();
                    throw new Exception("Not DataProducerException");
                }
            } else {
                if (dp.getReadCount() <= 0) {
                    throw new Exception("Empty");
                }
            }
        }
    }
}
