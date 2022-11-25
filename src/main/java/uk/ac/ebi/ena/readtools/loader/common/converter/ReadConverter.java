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
package uk.ac.ebi.ena.readtools.loader.common.converter;

import java.io.IOException;
import java.io.InputStream;

import uk.ac.ebi.ena.readtools.common.reads.QualityNormalizer;
import uk.ac.ebi.ena.readtools.loader.common.writer.ReadWriter;
import uk.ac.ebi.ena.readtools.loader.fastq.Read;

/**
 * Similar to {@link AutoNormalizeQualityReadConverter}, but here, base quality normalizer is provided explicitly. If no normalizer
 * is provided then normalization is not performed.
 */
public class ReadConverter extends AbstractReadConverter<Read> {

    private final ReadReader readReader;

    public ReadConverter(InputStream istream, ReadWriter<Read, ?> writer, String defaultReadIndex) {
        super(istream, writer);

        readReader = new ReadReader(defaultReadIndex);
    }

    public ReadConverter(InputStream istream, ReadWriter<Read, ?> writer, Long readLimit, String defaultReadIndex) {
        super(istream, writer, readLimit);

        readReader = new ReadReader(defaultReadIndex);
    }

    public ReadConverter(
            InputStream istream, ReadWriter<Read, ?> writer, QualityNormalizer normalizer, String defaultReadIndex) {

        super(istream, writer);

        readReader = new ReadReader(normalizer, defaultReadIndex);
    }

    public ReadConverter(
            InputStream istream,
            ReadWriter<Read, ?> writer,
            Long readLimit,
            QualityNormalizer normalizer,
            String defaultReadIndex) {

        super(istream, writer, readLimit);

        readReader = new ReadReader(normalizer, defaultReadIndex);
    }

    @Override
    public Read convert(InputStream inputStream) throws IOException {
        return readReader.read(inputStream);
    }
}
