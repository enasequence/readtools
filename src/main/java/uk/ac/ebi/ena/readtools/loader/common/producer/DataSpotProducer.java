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
package uk.ac.ebi.ena.readtools.loader.common.producer;

import java.io.IOException;
import java.io.InputStream;

import uk.ac.ebi.ena.readtools.common.reads.QualityNormalizer;
import uk.ac.ebi.ena.readtools.loader.fastq.DataSpot;

public class DataSpotProducer extends AbstractDataProducer<DataSpot> {

    private final DataSpotReader dataSpotReader;

    public DataSpotProducer(InputStream istream, QualityNormalizer normalizer, String defaultAttr) {
        super(istream);

        dataSpotReader = new DataSpotReader(normalizer, defaultAttr);
    }

    @Override
    public DataSpot produce(InputStream inputStream) throws IOException {
        return dataSpotReader.read(inputStream);
    }
}
