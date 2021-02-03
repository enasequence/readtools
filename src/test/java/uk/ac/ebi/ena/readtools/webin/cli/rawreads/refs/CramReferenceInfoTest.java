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
package uk.ac.ebi.ena.readtools.webin.cli.rawreads.refs;

import java.io.File;
import java.io.IOException;
import java.net.URL;
import java.net.URLDecoder;
import java.util.Map;

import org.junit.Assert;
import org.junit.Ignore;
import org.junit.Test;

import uk.ac.ebi.ena.readtools.webin.cli.rawreads.BamScannerTest;

public class 
CramReferenceInfoTest 
{
    @Ignore("Tests fails with an invalid sequence name error. The cram file that is being tested must be replaced with another one that has valid sequence names.")
    @Test public void
    testUnalignedCram() throws IOException, NoSuchFieldException, IllegalAccessException {
        URL url = BamScannerTest.class.getClassLoader().getResource( "uk/ac/ebi/ena/webin/cli/rawreads/nc-RNAs_DKC1_WT_1_1st_read.cram" );
        File file = new File( URLDecoder.decode( url.getFile(), "UTF-8" ) );
        CramReferenceInfo cri = new CramReferenceInfo();
        Map<?,?> result = cri.confirmFileReferences( file );
        Assert.assertEquals( 0, result.size() );
    }
    
    
    @Test public void
    testCorrectCramHeader() throws IOException
    {
        URL url = BamScannerTest.class.getClassLoader().getResource( "uk/ac/ebi/ena/webin/cli/rawreads/15194_1#135.cram" );
        File file = new File( URLDecoder.decode( url.getFile(), "UTF-8" ) );
        CramReferenceInfo cri = new CramReferenceInfo();
        Map<?,?> result = cri.confirmFileReferences( file );
        Assert.assertEquals( 66, result.size() );
    }
}
