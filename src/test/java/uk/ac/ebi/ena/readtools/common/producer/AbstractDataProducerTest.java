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

import java.io.EOFException;
import java.io.IOException;
import java.io.InputStream;
import java.time.Duration;
import java.time.LocalDateTime;

import org.junit.Assert;
import org.junit.Test;

import uk.ac.ebi.ena.readtools.loader.common.consumer.DataConsumer;
import uk.ac.ebi.ena.readtools.loader.common.consumer.DataConsumerException;
import uk.ac.ebi.ena.readtools.loader.common.consumer.Spot;
import uk.ac.ebi.ena.readtools.loader.common.producer.AbstractDataProducer;

public class AbstractDataProducerTest {

    @Test
    public void testRunDuration() throws InterruptedException {
        long waitTimeBeforeDummySpotMillis = 1;
        int dummySpotProduceCount = 5000;
        long expectedRunDurationSec = 1;

        AbstractDataProducer adp = new AbstractDataProducer(null, Duration.ofSeconds(expectedRunDurationSec)) {
            @Override
            public Spot produce(InputStream inputStream) throws IOException {
                //To prevent the test from running indefinitely (just in case).
                if (getReadCount() > dummySpotProduceCount) {
                    throw new EOFException();
                }

                try {
                    Thread.sleep(waitTimeBeforeDummySpotMillis);
                } catch (InterruptedException e) {
                    throw new RuntimeException(e);
                }

                return new Spot() {
                    @Override
                    public String getName() {
                        return "dummy";
                    }

                    @Override
                    public long getBaseCount() {
                        return 99;
                    }
                };
            }
        };
        adp.setConsumer(new DataConsumer() {
            @Override
            public void cascadeErrors() throws DataConsumerException {

            }

            @Override
            public void consume(Spot spot) throws DataConsumerException {

            }

            @Override
            public void setConsumer(DataConsumer dataConsumer) {

            }

            @Override
            public boolean isOk() {
                return true;
            }
        });

        LocalDateTime before = LocalDateTime.now();

        adp.start();
        adp.join();

        LocalDateTime after = LocalDateTime.now();

        long actualRunDurationSec = Duration.between(before, after).getSeconds();

        Assert.assertTrue(actualRunDurationSec >= expectedRunDurationSec);
        Assert.assertTrue(adp.getReadCount() > 0 && adp.getReadCount() < dummySpotProduceCount);
    }
}
