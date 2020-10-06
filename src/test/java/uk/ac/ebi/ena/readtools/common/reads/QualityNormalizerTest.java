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
package uk.ac.ebi.ena.readtools.common.reads;

import org.junit.Assert;
import org.junit.Test;
import uk.ac.ebi.ena.readtools.common.reads.normalizers.htsjdk.IlluminaQualityNormalizer;
import uk.ac.ebi.ena.readtools.common.reads.normalizers.htsjdk.SolexaQualityNormalizer;
import uk.ac.ebi.ena.readtools.common.reads.normalizers.htsjdk.StandardQualityNormalizer;

import java.nio.charset.StandardCharsets;

public class QualityNormalizerTest {

    @Test
    public void testStandardQualityNormalizer() {

        String quals = "!\"}~";
        byte[] expectedPhred = new byte[] {0, 1, 92, 93};

        byte[] actualPhred = quals.getBytes(StandardCharsets.UTF_8);

        new StandardQualityNormalizer().normalize(actualPhred);

        Assert.assertArrayEquals(expectedPhred, actualPhred);
    }

    @Test
    public void testSolexaQualityNormalizer() {

        String quals = ";<=>?@ABC{|}~";
        byte[] expectedPhred = new byte[] {0, 0, 0, 0, 0, 3, 4, 4, 5, 59, 60, 61, 62};

        byte[] actualPhred = quals.getBytes(StandardCharsets.UTF_8);

        new SolexaQualityNormalizer().normalize(actualPhred);

        Assert.assertArrayEquals(expectedPhred, actualPhred);
    }

    @Test
    public void testIlluminaQualityNormalizer() {

        String quals = "@A}~";
        byte[] expectedPhred = new byte[] {0, 1, 61, 62};

        byte[] actualPhred = quals.getBytes(StandardCharsets.UTF_8);

        new IlluminaQualityNormalizer().normalize(actualPhred);

        Assert.assertArrayEquals(expectedPhred, actualPhred);
    }
}
