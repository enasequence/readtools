package uk.ac.ebi.ena.readtools.old.webin.cli.rawreads;

import java.io.File;
import java.io.IOException;
import java.net.URL;
import java.net.URLDecoder;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.concurrent.atomic.AtomicBoolean;

import org.junit.Assert;
import org.junit.Test;

import uk.ac.ebi.ena.webin.cli.validator.message.ValidationResult;


public class 
BamScannerTest
{
    static class 
    MyScanner extends BamScanner
    { 
        protected void 
        logProcessedReadNumber( long cnt ) 
        {
            System.out.println( cnt );
        }
    };

    
    @Test public void
    testCorrectCram() throws IOException
    {
        URL url = BamScannerTest.class.getClassLoader().getResource( "uk/ac/ebi/ena/webin/cli/rawreads/18045_1#93.cram" );
        File file = new File( URLDecoder.decode( url.getFile(), "UTF-8" ) );
        
        BamScanner bs = new MyScanner();
        AtomicBoolean paired = new AtomicBoolean();
        RawReadsFile rf = new RawReadsFile();
        rf.setFilename( file.getPath() );
        ValidationResult vr = new ValidationResult();
        bs.readCramFile( vr, rf, paired );
        Assert.assertTrue( vr.isValid() );
        Assert.assertTrue( "cram should be paired", paired.get() );
    }
    
    
    @Test public void
    testIncorrectCram() throws IOException
    {
        URL url = BamScannerTest.class.getClassLoader().getResource( "uk/ac/ebi/ena/webin/cli/rawreads/15194_1#135.cram" );
        File file = new File( URLDecoder.decode( url.getFile(), "UTF-8" ) );
        BamScanner bs = new MyScanner();
        AtomicBoolean paired = new AtomicBoolean();
        RawReadsFile rf = new RawReadsFile();
        rf.setFilename( file.getPath() );
        ValidationResult vr = new ValidationResult();
        bs.readCramFile( vr, rf, paired );
        Assert.assertFalse( vr.isValid() );
    }

    
    @Test public void
    testIncorrectBAM() throws IOException
    {
        URL url = BamScannerTest.class.getClassLoader().getResource( "uk/ac/ebi/ena/webin/cli/rawreads/m54097_170904_165950.subreads.bam" );
        File file = new File( URLDecoder.decode( url.getFile(), "UTF-8" ) );
        BamScanner bs = new MyScanner();
        AtomicBoolean paired = new AtomicBoolean();
        RawReadsFile rf = new RawReadsFile();
        rf.setFilename( file.getPath() );
        ValidationResult vr = new ValidationResult();
        bs.readBamFile( vr, rf, paired );
        Assert.assertFalse( vr.isValid() );
    }
    
    
    
    @Test public void
    testCorrectBAM() throws IOException
    {
        URL url = BamScannerTest.class.getClassLoader().getResource( "uk/ac/ebi/ena/webin/cli/rawreads/OUTO500m_MetOH_narG_OTU18.bam" );
        Path file = Paths.get( new File( url.getFile() ).getCanonicalPath() );
        BamScanner bs = new MyScanner();
        AtomicBoolean paired = new AtomicBoolean();
        RawReadsFile rf = new RawReadsFile();
        rf.setFilename( String.valueOf( file ) );
        ValidationResult vr = new ValidationResult();
        bs.readBamFile( vr, rf, paired );
        Assert.assertTrue( vr.isValid() );
        Assert.assertFalse( paired.get() );
    }
    
    
    @Test public void
    testIncorrectBAMHeader() throws IOException
    {
        URL url = BamScannerTest.class.getClassLoader().getResource( "uk/ac/ebi/ena/webin/cli/rawreads/invalid.bam" );
        File file = new File( URLDecoder.decode( url.getFile(), "UTF-8" ) );
        BamScanner bs = new MyScanner();
        AtomicBoolean paired = new AtomicBoolean();
        RawReadsFile rf = new RawReadsFile();
        rf.setFilename( file.getPath() );
        Path ofile = Files.createTempFile( "BAMSCANNENR", "TEST" );
        ValidationResult vr = new ValidationResult( ofile.toFile() );
        bs.readCramFile( vr, rf, paired );
        Assert.assertFalse( vr.isValid() );
        Assert.assertTrue( new String( Files.readAllBytes( ofile ), StandardCharsets.UTF_8 ), new String( Files.readAllBytes( ofile ), StandardCharsets.UTF_8 ).contains( "File contains no valid reads" ) );
        
    }
}
