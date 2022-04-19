package uk.ac.ebi.ena.readtools.common.producer;

import org.junit.Test;
import uk.ac.ebi.ena.readtools.loader.common.consumer.DataConsumer;
import uk.ac.ebi.ena.readtools.loader.common.consumer.DataConsumerException;
import uk.ac.ebi.ena.readtools.loader.common.consumer.Spot;
import uk.ac.ebi.ena.readtools.loader.common.producer.AutoNormalizerDataSpotProducer;
import uk.ac.ebi.ena.readtools.loader.common.producer.DataProducerException;
import uk.ac.ebi.ena.readtools.loader.fastq.DataSpot;
import uk.ac.ebi.ena.readtools.utils.Utils;

import java.nio.file.Path;
import java.nio.file.Paths;
import java.time.Duration;

public class AutoNormalizerDataSpotProducerTest {
    private static final long FASTQ_VALIDATION_MAX_DURATION_MS = 4_000;

    @Test (timeout = FASTQ_VALIDATION_MAX_DURATION_MS + 1_000)
    public void test() throws Exception {
//        Path filePath = Paths.get("src/test/resources/tst.fastq.bz2");
//        Path filePath = Paths.get("src/test/resources/ERS10423544.fastq.gz");
        Path filePath = Paths.get("src/test/resources/ERS10423544.fastq.bz2");

        try (Utils.InputStreamFuture isf = Utils.openFastqInputStream(filePath)) {
            Duration duration = Duration.ofMillis(FASTQ_VALIDATION_MAX_DURATION_MS);
            AutoNormalizerDataSpotProducer dp =
                    new AutoNormalizerDataSpotProducer(isf.inputStream, duration, "", filePath.toString());
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
