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

import static org.junit.Assert.*;

import java.io.*;
import java.net.URL;
import java.nio.charset.StandardCharsets;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;

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

        InputStream bz2res = Utils.openFastqInputStream(filePath);
        BufferedReader br = new BufferedReader(new InputStreamReader(bz2res, StandardCharsets.UTF_8));

        List<String> sl = new ArrayList<>();

        while (true) {
            String line = br.readLine();
            if (null == line) {
                break;
            }
            sl.add(line);
        }
        bz2res.close();

        System.out.print(sl);

        assertEquals(32, sl.size());
    }

    @Test
    public void bz2StreamTryResource() {
        Path filePath = Paths.get("src/test/resources/tst.fastq.bz2");

        try (InputStream bz2res = Utils.openFastqInputStream(filePath)) {
            BufferedReader br = new BufferedReader(new InputStreamReader(bz2res, StandardCharsets.UTF_8));

            System.out.println(br.readLine());
        } catch (Exception exception) {
            exception.printStackTrace();
        }

        System.out.println("\n=====================================");

        try (InputStream bz2res = Utils.openFastqInputStream(filePath)) {
            BufferedReader br = new BufferedReader(new InputStreamReader(bz2res, StandardCharsets.UTF_8));
            List<String> sl = new ArrayList<>();

            while (true) {
                String line = br.readLine();
                if (null == line) {
                    break;
                }
                sl.add(line);
            }
            System.out.print(sl);

            assertEquals(32, sl.size());
        } catch (Exception exception) {
            exception.printStackTrace();
        }
    }
}
