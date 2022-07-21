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

import htsjdk.samtools.util.FastqQualityFormat;
import uk.ac.ebi.ena.readtools.common.reads.QualityNormalizer;
import uk.ac.ebi.ena.readtools.loader.fastq.Read;
import uk.ac.ebi.ena.readtools.utils.Utils;

import java.io.IOException;
import java.io.InputStream;

/**
 * Similar to {@link ReadConverter}. except here, quality normalization strategy is chosen automatically.
 */
public class FastqReadReadConverter extends AbstractReadConverter<Read> {

    private final String filePath;

    private final String defaultAttr;

    private volatile ReadReader readReader;

    public FastqReadReadConverter(InputStream istream, String defaultAttr, String filePath) {
        super(istream);

        this.defaultAttr = defaultAttr;
        this.filePath = filePath;
    }

    /**
     *
     * @param istream
     * @param readLimit Only read limited amount of reads.
     * @param defaultAttr
     * @param filePath
     */
    public FastqReadReadConverter(InputStream istream, Long readLimit, String defaultAttr, String filePath) {
        super(istream, readLimit);

        this.defaultAttr = defaultAttr;
        this.filePath = filePath;
    }

    @Override
    protected void begin() {
        try {
            FastqQualityFormat qualityFormat = Utils.detectFastqQualityFormat(filePath, null);

            QualityNormalizer normalizer = Utils.getQualityNormalizer(qualityFormat);

            readReader = new ReadReader(normalizer, defaultAttr);
        } catch (Exception ex) {
            throw new ConverterException(ex);
        }
    }

    @Override
    public Read convert(InputStream inputStream) throws IOException {
        return readReader.read(inputStream);
    }
}
