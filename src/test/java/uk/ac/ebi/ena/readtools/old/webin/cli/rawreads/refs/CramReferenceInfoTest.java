package uk.ac.ebi.ena.readtools.old.webin.cli.rawreads.refs;

import java.io.File;
import java.io.IOException;
import java.net.URL;
import java.net.URLDecoder;
import java.util.Map;

import org.junit.Assert;
import org.junit.Ignore;
import org.junit.Test;

import uk.ac.ebi.ena.readtools.old.webin.cli.rawreads.BamScannerTest;

public class 
CramReferenceInfoTest 
{
    @Ignore("Fails with an invalid sequence name exception. Sequence name validation rules have been updated in updated versions of htsjdk.")
    @Test public void
    testUnalignedCram() throws IOException
    {
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
