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

import java.io.*;
import java.net.URL;
import java.nio.charset.StandardCharsets;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import org.junit.Assert;
import org.junit.Test;

import htsjdk.samtools.fastq.FastqReader;

public class UtilsTest {

    @Test
    public void testSingleFastqTest() throws IOException {
        String expectedUpper = "GACGGATCTATAGCAAAACT";
        String expectedLower = "gacggatctatagcaaaact";

        File output = File.createTempFile("FASTQ", "FASTQ");
        output.delete();

        URL url = UtilsTest.class.getClassLoader().getResource("uracil-bases.fastq");

        Utils.replaceUracilBasesInFastq(url.getFile(), output.getAbsolutePath());

        FastqReader reader = new FastqReader(output);

        String recReadString1 = reader.next().getReadString();
        String recReadString2 = reader.next().getReadString();

        reader.close();

        Assert.assertEquals(expectedUpper, recReadString1);
        Assert.assertEquals(expectedLower, recReadString2);
    }

    @Test
    public void bz2Stream() throws IOException {
        Path filePath = Paths.get("src/test/resources/tst.fastq.bz2");
        Utils.InputStreamFuture bz2res = Utils.createBz2InputStream(filePath);
        BufferedReader br = new BufferedReader(new InputStreamReader(bz2res.inputStream, StandardCharsets.UTF_8));

        System.out.println(br.readLine());

        Integer exitCode = Utils.getBz2ReadResult(bz2res.future);
        System.out.println("exitCode " + exitCode);

        bz2res = Utils.createBz2InputStream(filePath);
        br = new BufferedReader(new InputStreamReader(bz2res.inputStream, StandardCharsets.UTF_8));
        List<String> sl = new ArrayList<>();

        System.out.println("\n=====================================");

        while (true) {
            String line = br.readLine();
            if (null == line) {
                break;
            }
            sl.add(line);
        }

        exitCode = Utils.getBz2ReadResult(bz2res.future);
        System.out.println("exitCode " + exitCode);
        System.out.print(sl);
    }
}
