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
package uk.ac.ebi.ena.readtools.webin.cli.rawreads;

import java.io.File;
import java.io.IOException;

import org.junit.Test;

public class SamScannerTest {
    @Test public void
    test() throws IOException {
        File file = new File("/Users/eugene/temp/samtest/mixed_reads.cram");
        File refFile = new File("/Users/eugene/temp/samtest/dummy_ref.fa");

        SamScanner samScanner = new SamScanner();
//        samScanner.checkSamFile(file, refFile);
        samScanner.checkSamFile2(file, refFile);
    }
}
