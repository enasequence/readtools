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
package uk.ac.ebi.ena.readtools.common.converter;

import java.io.EOFException;
import java.io.IOException;
import java.io.InputStream;

import org.junit.Assert;
import org.junit.Test;

import uk.ac.ebi.ena.readtools.loader.common.converter.AbstractReadConverter;
import uk.ac.ebi.ena.readtools.loader.common.writer.ReadWriter;
import uk.ac.ebi.ena.readtools.loader.common.writer.ReadWriterException;
import uk.ac.ebi.ena.readtools.loader.common.writer.Spot;

public class AbstractReadConverterTest {

    @Test
    public void testRunLimit() {
        long readLimit = 10;

        AbstractReadConverter converter = new AbstractReadConverter(
                null,
                new ReadWriter() {
                    @Override
                    public void cascadeErrors() throws ReadWriterException {

                    }

                    @Override
                    public void write(Spot spot) throws ReadWriterException {

                    }

                    @Override
                    public void setWriter(ReadWriter readWriter) {

                    }
                },
                readLimit) {
            @Override
            public Spot convert(InputStream inputStream) throws IOException {
                //To prevent the test from running indefinitely (just in case).
                if (getReadCount() > readLimit + 1) {
                    throw new EOFException();
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

                    @Override
                    public long getSizeBytes() {
                        return 1024 * 99L;
                    }
                };
            }
        };

        converter.run();

        Assert.assertEquals(readLimit, converter.getReadCount());
    }
}
