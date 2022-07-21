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

import uk.ac.ebi.ena.readtools.common.reads.QualityNormalizer;
import uk.ac.ebi.ena.readtools.loader.fastq.Read;

import java.io.IOException;
import java.io.InputStream;

/**
 * Similar to {@link FastqReadReadConverter}, but here, base quality normalizer is provided explicitly. If no normalizer
 * is provided then normalization is not performed.
 */
public class ReadConverter extends AbstractReadConverter<Read> {

    private final ReadReader readReader;

    public ReadConverter(InputStream istream, String defaultAttr) {
        super(istream);

        readReader = new ReadReader(defaultAttr);
    }

    public ReadConverter(InputStream istream, Long readLimit, String defaultAttr) {
        super(istream, readLimit);

        readReader = new ReadReader(defaultAttr);
    }

    public ReadConverter(InputStream istream, QualityNormalizer normalizer, String defaultAttr) {
        super(istream);

        readReader = new ReadReader(normalizer, defaultAttr);
    }

    public ReadConverter(InputStream istream, Long readLimit, QualityNormalizer normalizer, String defaultAttr) {
        super(istream, readLimit);

        readReader = new ReadReader(normalizer, defaultAttr);
    }

    @Override
    public Read convert(InputStream inputStream) throws IOException {
        return readReader.read(inputStream);
    }
}
