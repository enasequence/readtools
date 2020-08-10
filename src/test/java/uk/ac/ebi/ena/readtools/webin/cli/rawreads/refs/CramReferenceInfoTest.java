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
package uk.ac.ebi.ena.readtools.webin.cli.rawreads.refs;

import java.io.File;
import java.io.IOException;
import java.lang.reflect.Field;
import java.lang.reflect.Modifier;
import java.net.URL;
import java.net.URLDecoder;
import java.util.Map;
import java.util.regex.Pattern;

import org.junit.Assert;
import org.junit.Test;

import htsjdk.samtools.SAMSequenceRecord;

import uk.ac.ebi.ena.readtools.webin.cli.rawreads.BamScannerTest;

public class 
CramReferenceInfoTest 
{
    @Test public void
    testUnalignedCram() throws IOException, NoSuchFieldException, IllegalAccessException {
        //Replace built-in sequence name validation pattern with a dummy one as it's not possible to skip/turn off validation.
        Field field = SAMSequenceRecord.class.getDeclaredField("LEGAL_RNAME_PATTERN");
        field.setAccessible(true);

        Field modifiersField = Field.class.getDeclaredField("modifiers");
        modifiersField.setAccessible(true);
        modifiersField.setInt(field, field.getModifiers() & ~Modifier.FINAL);

        Object oldPattern = field.get(null);
        field.set(null, Pattern.compile("."));

        URL url = BamScannerTest.class.getClassLoader().getResource( "uk/ac/ebi/ena/webin/cli/rawreads/nc-RNAs_DKC1_WT_1_1st_read.cram" );
        File file = new File( URLDecoder.decode( url.getFile(), "UTF-8" ) );
        CramReferenceInfo cri = new CramReferenceInfo();
        Map<?,?> result = cri.confirmFileReferences( file );
        Assert.assertEquals( 0, result.size() );

        //Set the old pattern back in.
        field.set(null, oldPattern);
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
