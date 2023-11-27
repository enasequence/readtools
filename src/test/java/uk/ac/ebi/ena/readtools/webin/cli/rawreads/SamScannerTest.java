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
import java.io.UnsupportedEncodingException;
import java.net.URL;
import java.net.URLDecoder;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.time.Duration;
import java.time.LocalDateTime;
import java.util.concurrent.atomic.AtomicBoolean;

import org.junit.Assert;
import org.junit.Test;

import uk.ac.ebi.ena.readtools.refactored.validator.ReadsValidationException;
import uk.ac.ebi.ena.readtools.refactored.validator.ValidatorWrapper;
import uk.ac.ebi.ena.webin.cli.validator.message.ValidationResult;


public class
SamScannerTest
{
    private static final Long READ_LIMIT = 1024L;

    static class 
    MyScanner extends SamScanner
    {
        public MyScanner() {
            super(READ_LIMIT);
        }

        protected void
        logProcessedReadNumber( long cnt ) 
        {
            System.out.println( cnt );
        }
    };



    
    @Test public void
    testCorrectCram() throws IOException
    {
        URL url = SamScannerTest.class.getClassLoader().getResource("rawreads/18045_1#93.cram");
        File file = new File( URLDecoder.decode( url.getFile(), "UTF-8" ) );
        
        SamScanner bs = new MyScanner();
        AtomicBoolean paired = new AtomicBoolean();
        RawReadsFile rf = new RawReadsFile();
        rf.setFilename( file.getPath() );
        ValidationResult vr = new ValidationResult();
        bs.readCramFile( vr, rf, paired, true );
        Assert.assertTrue( vr.isValid() );
        Assert.assertTrue( "cram should be paired", paired.get() );
    }
    
    
    @Test public void
    testIncorrectCram() throws IOException
    {
        URL url = SamScannerTest.class.getClassLoader().getResource("rawreads/15194_1#135.cram");
        File file = new File( URLDecoder.decode( url.getFile(), "UTF-8" ) );
        SamScanner bs = new MyScanner();
        AtomicBoolean paired = new AtomicBoolean();
        RawReadsFile rf = new RawReadsFile();
        rf.setFilename( file.getPath() );
        ValidationResult vr = new ValidationResult();
        bs.readCramFile( vr, rf, paired, true );
        Assert.assertFalse( vr.isValid() );
    }

    
    @Test public void
    testIncorrectBAM() throws IOException
    {
        URL url = SamScannerTest.class.getClassLoader().getResource("rawreads/m54097_170904_165950.subreads.bam");
        File file = new File( URLDecoder.decode( url.getFile(), "UTF-8" ) );
        SamScanner bs = new MyScanner();
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
        URL url = SamScannerTest.class.getClassLoader().getResource("rawreads/OUTO500m_MetOH_narG_OTU18.bam");
        Path file = Paths.get( new File( url.getFile() ).getCanonicalPath() );
        SamScanner bs = new MyScanner();
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
        URL url = SamScannerTest.class.getClassLoader().getResource("rawreads/invalid.bam");
        File file = new File( URLDecoder.decode( url.getFile(), "UTF-8" ) );
        SamScanner bs = new MyScanner();
        AtomicBoolean paired = new AtomicBoolean();
        RawReadsFile rf = new RawReadsFile();
        rf.setFilename( file.getPath() );
        Path ofile = Files.createTempFile( "BAMSCANNENR", "TEST" );
        ValidationResult vr = new ValidationResult( ofile.toFile() );
        bs.readCramFile( vr, rf, paired, true );
        Assert.assertFalse( vr.isValid() );
        Assert.assertTrue( new String( Files.readAllBytes( ofile ), StandardCharsets.UTF_8 ), new String( Files.readAllBytes( ofile ), StandardCharsets.UTF_8 ).contains( "File contains no valid reads" ) );
        
    }

    @Test
    public void testRunDuration() throws IOException {
        URL url = SamScannerTest.class.getClassLoader().getResource("bam2fastq/3fastq/M2241_BLV_sense.bam");
        Path file = Paths.get( new File( url.getFile() ).getCanonicalPath() );

        //Adjust this if the file above is changed.
        long fileReadCount = 129;

        //Lower further if scanner completes faster than this.
        Duration expectedRunDuration = Duration.ofMillis(1);

        RawReadsFile rf = new RawReadsFile();
        rf.setFilename( String.valueOf( file ) );

        AtomicBoolean paired = new AtomicBoolean();

        ValidationResult vr = new ValidationResult();

        long processedReadCount[] = {0};
        SamScanner bs = new SamScanner(READ_LIMIT, expectedRunDuration) {
            @Override
            protected void logProcessedReadNumber(long cnt) {
                processedReadCount[0] = cnt;
            }
        };

        LocalDateTime before = LocalDateTime.now();

        bs.readBamFile( vr, rf, paired );

        LocalDateTime after = LocalDateTime.now();

        Duration actualRunDuration = Duration.between(before, after);

        Assert.assertTrue(actualRunDuration.getNano() >= expectedRunDuration.getNano());

        //If the scanner returned within the the time budget as it was supposed to
        //then processed read count should be lower than total read count in the file.
        //Which means the timeout functionality worked.
        Assert.assertTrue(processedReadCount[0] > 0 && processedReadCount[0] < fileReadCount);
    }

    @Test
    public void srv1() throws UnsupportedEncodingException, ReadsValidationException {
        URL url = SamScannerTest.class.getClassLoader().getResource("rawreads/15194_1#135.cram");
        File file = new File(URLDecoder.decode(url.getFile(), "UTF-8"));
        ValidatorWrapper.validateSam(file);
    }
    @Test
    public void srv2() throws UnsupportedEncodingException, ReadsValidationException {
        URL url = SamScannerTest.class.getClassLoader().getResource("rawreads/m54097_170904_165950.subreads.bam");
        File file = new File( URLDecoder.decode(url.getFile(), "UTF-8"));
        ValidatorWrapper.validateSam(file);
    }
}
