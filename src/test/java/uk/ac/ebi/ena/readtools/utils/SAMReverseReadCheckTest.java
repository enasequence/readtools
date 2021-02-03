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
package uk.ac.ebi.ena.readtools.utils;

import java.io.File;
import java.net.URL;

import org.junit.Assert;
import org.junit.Test;

import uk.ac.ebi.ena.readtools.sam.Sam2FastqTest;

public class SAMReverseReadCheckTest {

    @Test
    public void reverseReadPresent() {
        URL url = Sam2FastqTest.class.getClassLoader().getResource("sam2fastq/reverse-read.sam");
        File file = new File( url.getFile() );

        SAMReverseReadCheck.Output out = new SAMReverseReadCheck(file).check();

        Assert.assertEquals("rid-0001", out.readName);
        Assert.assertEquals("TCGATCGATCGATCGA", out.bases);
        Assert.assertEquals("PONMLKJIHGFEDCBA", out.qualities);
        Assert.assertEquals(157, out.flags);
    }

    @Test
    public void reverseReadNotPresent() {
        URL url = Sam2FastqTest.class.getClassLoader().getResource("bam2fastq/2fastq/8855_124.bam");
        File file = new File( url.getFile() );

        Assert.assertNull(new SAMReverseReadCheck(file).check());
    }
}
