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
package uk.ac.ebi.ena.readtools.fastq;

import org.junit.Assert;
import org.junit.Test;

import java.io.File;
import java.io.IOException;
import java.net.URISyntaxException;
import java.net.URL;
import java.nio.file.Files;
import java.nio.file.Paths;

public class Fastq2SamTest {

    @Test
    public void testPairedFastq() throws IOException, URISyntaxException {
        File output = File.createTempFile( "SAM", ".sam" );
        output.delete();

        URL url1 = Fastq2SamTest.class.getClassLoader().getResource( "fastq2sam/fastq2sam1.fastq");
        URL url2 = Fastq2SamTest.class.getClassLoader().getResource( "fastq2sam/fastq2sam2.fastq");

        Fastq2Sam fastq2Sam = new Fastq2Sam("sample001", url1.getFile(), url2.getFile(), output.getAbsolutePath(), null);
        fastq2Sam.convert();

        byte[] expected = Files.readAllBytes(Paths.get(Fastq2SamTest.class.getClassLoader().getResource("fastq2sam/fastq2sam.sam").toURI()));
        byte[] actual = Files.readAllBytes(output.toPath());

        Assert.assertArrayEquals(expected, actual);
    }
}
