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
import java.time.Duration;

import org.junit.Test;

import uk.ac.ebi.ena.readtools.loader.common.consumer.DataConsumer;
import uk.ac.ebi.ena.readtools.loader.common.consumer.DataConsumerException;
import uk.ac.ebi.ena.readtools.loader.common.consumer.Spot;
import uk.ac.ebi.ena.readtools.loader.common.producer.AutoNormalizerDataSpotProducer;
import uk.ac.ebi.ena.readtools.loader.common.producer.DataProducerException;
import uk.ac.ebi.ena.readtools.loader.fastq.DataSpot;
import uk.ac.ebi.ena.readtools.utils.Utils;

public class AutoNormalizerDataSpotProducerTest {
    private static final long FASTQ_VALIDATION_MAX_DURATION_MS = 400_000;

    @Test (timeout = FASTQ_VALIDATION_MAX_DURATION_MS + 1_000)
    public void test() throws Exception {
//        Path filePath = Paths.get("src/test/resources/tst.fastq.bz2");
//        Path filePath = Paths.get("src/test/resources/ERS10423544.fastq.gz");
        Path filePath = Paths.get("src/test/resources/ERS10423544.fastq.bz2");

        try (InputStream is = Utils.openFastqInputStream(filePath)) {
            Duration duration = Duration.ofMillis(FASTQ_VALIDATION_MAX_DURATION_MS);
            AutoNormalizerDataSpotProducer dp =
                    new AutoNormalizerDataSpotProducer(is, duration, "", filePath.toString());
            dp.setName(filePath.toFile().getName());
            dp.setConsumer(new DataConsumer<DataSpot, Spot>() {
                @Override
                public void cascadeErrors() throws DataConsumerException {
                }

                @Override
                public void
                consume(DataSpot spot) {
                }

                @Override
                public void setConsumer(DataConsumer<Spot, ?> dataConsumer) {
                    throw new RuntimeException("Not implemented");
                }

                @Override
                public boolean isOk() {
                    return true;
                }
            });
            dp.start();
            dp.join();
            if (!dp.isOk()) {
                if (dp.getStoredException() instanceof DataProducerException) {
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
