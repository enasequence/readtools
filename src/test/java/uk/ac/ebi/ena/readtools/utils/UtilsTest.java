/*
 * Copyright 2010-2020 EMBL - European Bioinformatics Institute
 * Licensed under the Apache License, Version 2.0 (the "License"); you may not use this
 * file except in compliance with the License. You may obtain a copy of the License at
 * http://www.apache.org/licenses/LICENSE-2.0
 * Unless required by applicable law or agreed to in writing, software distributed under the
 * License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
 * CONDITIONS OF ANY KIND, either express or implied. See the License for the
 * specific language governing permissions and limitations under the License.
 */
package uk.ac.ebi.ena.readtools.utils;

import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import org.junit.Assert;
import org.junit.Test;
import uk.ac.ebi.ena.readtools.fastq.Fastq2Sam;

import java.io.File;
import java.io.IOException;
import java.net.URL;

public class UtilsTest {

    @Test
    public void testSingleFastqTest() throws IOException {
        String expectedUpper = "GACGGATCTATAGCAAAACT";
        String expectedLower = "gacggatctatagcaaaact";

        File output = File.createTempFile( "FASTQ", "FASTQ" );
        output.delete();

        URL url = UtilsTest.class.getClassLoader().getResource( "uracil-bases.fastq");

        Utils.replaceUracilBasesInFastq(url.getFile(), output.getAbsolutePath());

        FastqReader reader = new FastqReader(output);

        String recReadString1 = reader.next().getReadString();
        String recReadString2 = reader.next().getReadString();

        reader.close();

        Assert.assertEquals(expectedUpper, recReadString1);
        Assert.assertEquals(expectedLower, recReadString2);
    }
}
