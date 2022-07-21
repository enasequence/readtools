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

import org.junit.Assert;
import org.junit.Test;
import uk.ac.ebi.ena.webin.cli.validator.message.ValidationMessage.Severity;
import uk.ac.ebi.ena.webin.cli.validator.message.ValidationResult;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.net.URL;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import java.util.concurrent.ThreadLocalRandom;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import java.util.zip.GZIPOutputStream;

public class
FastqScannerTest
{
    private static final int expected_reads = 10_000;
    
    static class
    MyScanner extends FastqScanner
    {
		public
		MyScanner( int expected_size )
		{
			super( expected_size );
		}

		
	    @Override protected void
	    logProcessedReadNumber( long count )
	    {
	        String msg = String.format( "\rProcessed %16d read(s)", count );
	        logFlushMsg( msg );
	    }

	    
	    @Override protected void
	    logFlushMsg( String msg )
	    {
	        System.out.print( msg );
	        System.out.flush();
	    }
    }

    
    @Test public void
    testSingle() throws Throwable
    {
        URL  url1 = FastqScannerTest.class.getClassLoader().getResource( "uk/ac/ebi/ena/webin/cli/rawreads/EP0_GTTCCTT_S1.txt.gz" );
        FastqScanner fs = new MyScanner( expected_reads );
        RawReadsFile rf = new RawReadsFile();
        
        rf.setFilename( new File( url1.getFile() ).getCanonicalPath() );

        ValidationResult vr = new ValidationResult();

        fs.checkFiles( vr, rf );

        Assert.assertTrue( vr.isValid() );
    }

    
    @Test public void
    testSingleDuplications() throws Throwable
    {
        URL  url1 = FastqScannerTest.class.getClassLoader().getResource( "uk/ac/ebi/ena/webin/cli/rawreads/EP0_GTTCCTT_S1.txt.dup.gz" );
        FastqScanner fs = new MyScanner( expected_reads );
        RawReadsFile rf = new RawReadsFile();
        
        rf.setFilename( new File( url1.getFile() ).getCanonicalPath() );

        ValidationResult vr = new ValidationResult();

        fs.checkFiles( vr, rf );
        
        Assert.assertEquals( 2, vr.count(Severity.ERROR) );
    }


    @Test public void
    testPaired() throws Throwable
    {
        URL  url1 = FastqScannerTest.class.getClassLoader().getResource( "uk/ac/ebi/ena/webin/cli/rawreads/EP0_GTTCCTT_0.txt.gz" );
        FastqScanner fs = new MyScanner( expected_reads );
        RawReadsFile rf = new RawReadsFile();
        
        rf.setFilename( new File( url1.getFile() ).getCanonicalPath() );

        ValidationResult vr = new ValidationResult();

        fs.checkFiles( vr, rf );
        
        Assert.assertTrue( vr.isValid() );
    }

    
    @Test public void
    testPair() throws Throwable
    {
        URL  url1 = FastqScannerTest.class.getClassLoader().getResource( "uk/ac/ebi/ena/webin/cli/rawreads/EP0_GTTCCTT_S1.txt.gz" );
        URL  url2 = FastqScannerTest.class.getClassLoader().getResource( "uk/ac/ebi/ena/webin/cli/rawreads/EP0_GTTCCTT_S2.txt.gz" );
        FastqScanner fs = new MyScanner( expected_reads );
        RawReadsFile rf1 = new RawReadsFile();
        rf1.setFilename( new File( url1.getFile() ).getCanonicalPath() );
        
        RawReadsFile rf2 = new RawReadsFile();
        rf2.setFilename( new File( url2.getFile() ).getCanonicalPath() );

        ValidationResult vr = new ValidationResult();

        fs.checkFiles( vr, rf1, rf2 );
        
        Assert.assertTrue( vr.isValid() );
    }


    
    private Path
    saveRandomized(String content, Path folder, boolean gzip, String... suffix) throws IOException
    {
        Path file = Files.createTempFile( "_content_", "_content_" );
        Files.write( file, content.getBytes( StandardCharsets.UTF_8 ), StandardOpenOption.TRUNCATE_EXISTING, StandardOpenOption.SYNC );
        Path path = Files.createTempFile( folder, "COPY", ( file.getName( file.getNameCount() - 1 ) + ( suffix.length > 0 ? Stream.of( suffix ).collect( Collectors.joining( ".", ".", "" ) ) : "" ) ));
        OutputStream os;
        Files.copy( file, ( os = gzip ? new GZIPOutputStream( Files.newOutputStream( path, StandardOpenOption.TRUNCATE_EXISTING, StandardOpenOption.SYNC ) ) 
                                      : Files.newOutputStream( path, StandardOpenOption.TRUNCATE_EXISTING, StandardOpenOption.SYNC ) ) );
        os.flush();
        os.close();
        Assert.assertTrue( Files.exists( path ) );
        Assert.assertTrue( Files.isRegularFile( path ) );
        return path;
    }


    private File
    createOutputFolder() throws IOException
    {
        File output = File.createTempFile( "test", "test" );
        Assert.assertTrue( output.delete() );
        Assert.assertTrue( output.mkdirs() );
        return output;
    }
    
    /* Test cases */
    /* 1. Single run with duplication */
    /* 2. Paired run with one file with duplication */
    /* 3. Paired run with two files with duplication in first file */
    /* 4. Paired run with two files with duplication from second file in second file */
    @Test public void 
    testCase1() throws Throwable
    {
        File output_dir = createOutputFolder();
        Path f1 = saveRandomized( "@NAME1/1\nACGT\n+\n1234\n"
                                + "@NAME1/1\nACGT\n+\n1234\n" 
                                + "@NAME2/1\nACGT\n+\n1234", output_dir.toPath(), true, "fastq-1", "gz" );
        FastqScanner fs = new MyScanner( expected_reads );
        RawReadsFile rf1 = new RawReadsFile();
        rf1.setFilename( f1.toFile().getCanonicalPath() );

        ValidationResult vr = new ValidationResult();

        fs.checkFiles( vr, rf1 );

        Assert.assertFalse( fs.getPaired() );
        Assert.assertEquals( 1, vr.count(Severity.ERROR) );
    }

    
    @Test public void 
    testCase2() throws Throwable
    {
        File output_dir = createOutputFolder();
        Path f1 = saveRandomized( "@NAME1/1\nACGT\n+\n1234\n"
                                + "@NAME1/2\nACGT\n+\n1234\n" 
                                + "@NAME1/1\nACGT\n+\n1234", output_dir.toPath(), true, "fastq-1", "gz" );
        FastqScanner fs = new MyScanner( expected_reads );
        RawReadsFile rf1 = new RawReadsFile();
        rf1.setFilename( f1.toFile().getCanonicalPath() );

        ValidationResult vr = new ValidationResult();

        fs.checkFiles( vr, rf1 );

        Assert.assertFalse( fs.getPaired() );
        Assert.assertEquals( 1, vr.count(Severity.ERROR) );
    }
    
    
    @Test public void 
    testCase3() throws Throwable
    {
        File output_dir = createOutputFolder();
        Path f1 = saveRandomized( "@NAME1\nACGT\n+\n1234\n"
                                + "@NAME1\nACGT\n+\n1234", output_dir.toPath(), true, "fastq-1", "gz" );
        Path f2 = saveRandomized( "@NAME2\nACGT\n+\n1234", output_dir.toPath(), true, "fastq-2", "gz" );
        FastqScanner fs = new MyScanner( expected_reads );
        RawReadsFile rf1 = new RawReadsFile();
        rf1.setFilename( f1.toFile().getCanonicalPath() );
        
        RawReadsFile rf2 = new RawReadsFile();
        rf2.setFilename( f2.toFile().getCanonicalPath() );

        ValidationResult vr = new ValidationResult();

        fs.checkFiles( vr, rf1, rf2 );
        
        Assert.assertEquals( 1, vr.count(Severity.ERROR) );
    }

    /* */
    @Test public void 
    testCase5() throws Throwable
    {
        File output_dir = createOutputFolder();
        Path f1 = saveRandomized( "@NAME1\nACGT\n+\n1234", output_dir.toPath(), true, "fastq-1", "gz" );
        Path f2 = saveRandomized( "@NAME2\nACGT\n+\n1234\n"
                                + "@NAME2\nACGT\n+\n2341", output_dir.toPath(), true, "fastq-2", "gz" );
        
        FastqScanner fs = new MyScanner( expected_reads );
        RawReadsFile rf1 = new RawReadsFile();
        rf1.setFilename( f1.toFile().getCanonicalPath() );
        
        RawReadsFile rf2 = new RawReadsFile();
        rf2.setFilename( f2.toFile().getCanonicalPath() );

        ValidationResult vr = new ValidationResult();

        fs.checkFiles( vr, rf1, rf2 );

        Assert.assertEquals( 1, vr.count(Severity.ERROR) );
    }

    
    /* three read labels */
    @Test public void 
    testCase6() throws Throwable
    {
        File output_dir = createOutputFolder();
        Path f1 = saveRandomized( "@NAME\nACGT\n+\n1234", output_dir.toPath(), true, "fastq-1", "gz" );
        Path f2 = saveRandomized( "@NAME/2\nACGT\n+\n1234\n"
                                + "@NAME/3\nACGT\n+\n2341", output_dir.toPath(), true, "fastq-2", "gz" );
        
        FastqScanner fs = new MyScanner( expected_reads );
        RawReadsFile rf1 = new RawReadsFile();
        rf1.setFilename( f1.toFile().getCanonicalPath() );
        
        RawReadsFile rf2 = new RawReadsFile();
        rf2.setFilename( f2.toFile().getCanonicalPath() );

        ValidationResult vr = new ValidationResult();

        fs.checkFiles( vr, rf1, rf2 );
        
        Assert.assertEquals( 1, vr.count(Severity.ERROR) );
    }


    /* Wrong pair set in two files */
    @Test public void 
    testCase7() throws Throwable
    {
        File output_dir = createOutputFolder();
        Path f1 = saveRandomized( "@NAME1/1\nACGT\n+\n1234", output_dir.toPath(), true, "fastq-1", "gz" );
        Path f2 = saveRandomized( "@NAME2/2\nACGT\n+\n1234\n"
                                + "@NAME/2\nACGT\n+\n2341", output_dir.toPath(), true, "fastq-2", "gz" );
        
        FastqScanner fs = new MyScanner( expected_reads );
        RawReadsFile rf1 = new RawReadsFile();
        rf1.setFilename( f1.toFile().getCanonicalPath() );
        
        RawReadsFile rf2 = new RawReadsFile();
        rf2.setFilename( f2.toFile().getCanonicalPath() );

        ValidationResult vr = new ValidationResult();

        fs.checkFiles( vr, rf1, rf2 );

        Assert.assertEquals( 1, vr.count(Severity.ERROR) );
    }
    
    
    /*  PacBio RS II Wrong pair set in one file */
    @Test public void 
    testCase9() throws Throwable
    {
        File output_dir = createOutputFolder();
        Path f1 = saveRandomized( 
        "@m150522_164945_42237_c100824932550000001823179511031571_s1_p0/134/0_15617\n"
      + "ACCAACAAGAGA\n"
      + "+\n"
      + "\"#'#&''$,-,)\n"
      
      + "@m150522_164945_42237_c100824932550000001823179511031571_s1_p0/135/0_15617\n"
      + "ACCAACAAGAGA\n"
      + "+\n"
      + "\"#'#&''$,-,)\n"
      
      + "@m150522_164945_42237_c100824932550000001823179511031571_s1_p0/134/15661_24597\n"
      + "ACCACCCTTAT\n"
      + "+\n"
      + "/&\"-('--.#/\n"

      + "@m150522_164945_42237_c100824932550000001823179511031571_s1_p0/135/15661_24597\n"
      + "ACCACCCTTAT\n"
      + "+\n"
      + "/&\"-('--.#/\n",
      	output_dir.toPath(), true, "fastq-9", "gz" );
        
        FastqScanner fs = new MyScanner( expected_reads );
        RawReadsFile rf1 = new RawReadsFile();
        rf1.setFilename( f1.toFile().getCanonicalPath() );

        ValidationResult vr = new ValidationResult();

        fs.checkFiles( vr, rf1 );

        Assert.assertFalse( fs.getPaired() );
        Assert.assertEquals( 0, vr.count(Severity.ERROR) );
    }
    
    
    /*  PacBio RS II Wrong pair set in one file */
    @Test public void 
    testCase10() throws Throwable
    {
        File output_dir = createOutputFolder();
        Path f1 = saveRandomized( 
            	"@CL100031435L2C001R001_63456 :TGGTCCTTCC\n"
        	  + "ATCCTTATAGTGGGCCAAGCTCTCACATGCAAACACGTTGCCTGCTGGATTTTGTTTCAGAGAGGAAATGTTTATGTGAGACAGAAAAAGCCGGGGGCC\n"
        	  + "+\n"
        	  + "FDFFFFFFFG3FFFFFF7?D8FFFFF=>DFG(FBEFFFFDF;FF?FF<8FGGFFFCBB8F0F@FBC?FAGFEE>.FFEFCF:?E(E@;3*(,FD/BFE-\n"
              + "@CL100031435L2C001R001_63472 :TGGTCCTTCA\n"
              + "AGCACTATGGAAAATGGAACCTCTTTGGGTCCTGTAATCAC\n"
              + "+\n"
              + "FFFFFFEGFFFGEFGFFFGFGGFEFDFFFFGGEFEFBG;FG\n"
              + "@CL100031435L2C001R001_63487 :TGGTCCTTCC\n"
              + "CCGCTCTGCCCACCACAATCCGGCCCCTGTGTACGGCAACACAGGGCCCAGGCAAGCAGATCCTTCCTGCTGGGAGCTCCAGCTTGTAGAATTTCACCC\n"
              + "+\n"
              + "FFFFFFFFFEFFFFFFD5FEFFEEEFFEBCE?EEFFFECD:DAFFCDFFBFFDB@>FCFBEEFF>F-CBFFF&E9F=BBF:4@B?B45E2F+A-EFFCA\n"
              + "@CL100031435L2C001R001_63499 :TGGTCCTTGC\n"
              + "CTCTTTTGATGCCTTCTGTATTGAGACACTTTCAGCAGAGTCCCAAGGCCAGGTGAAGAGGAAGCATGGCGGGGGCAGGTCGGGGGCCAACGGGTCTTAG\n"
              + "+\n"
              + "FEFF@E9F=FFF7FDD@FDFFFF4F;EFFFFAA8FEADDFFEFFE5F?FEBDEFFACD;BFFBFD2F>EDD:7FDDF'FDFF)+7FF1F1<4FF:@A@/F\n"
              + "@CL100031435L2C001R001_63541 :TGGTCCTTGC\n"
              + "CATGGGGGCATCTTCCTGCTCCAGGCACAAAGCTCTGGCTGTCACAACCCAG\n"
              + "+\n"
              + "FFFGFF>GFFGFFBFAGGFEFFDEFGFFEEAGGFFFBB>GGBFF6FFFGF@F\n"
              + "@CL100031435L2C001R001_63542 :TGGTCCTTGC\n"
              + "CCTCATTCTAAAATGCAGGTGCTCTGCTGACAGCAAAATTCTGTGTTTGAGCTGTGCAGTCTTAAACCAAGACGTGGGAGGCCAGGCGGTTTTAGAAAG\n"
              + "+\n"
              + "AAE>AFE4EFFA3FBAEDDFBEFDDA?F=;E?6?>FC4BF@EFC5A9FE=ECFCCD<CEEBF(EACE1F<5D(0(+DA;C9@1F&3,F76&B?DF1<30\n"
              + "@CL100031435L2C001R001_63548 :TGGTCCTTGC\n"
              + "GTGGTGGTGGTGTTTTTTTGGTTTGGGGTTTGGGTTTTTTTCAGGGAA\n"
              + "+\n"
              + "GFGFFFFFF=DGFFGGEFFFFFBGFF>FFFGFFFEGFGGEEGGGGFGF\n"
              + "@CL100031435L2C001R001_63553 :TGGTCCTTCC\n"
              + "GCACTTCTCGAGCTTCACATTCTAATGAGAACAATTTCCTTGGATTCATTGGTGTTGCCATTTTTTTGTTGACTCATTCAAAAAACAAATTAGCTGAGGT\n"
              + "+\n"
              + "FFFFFFFFFFDFEGFFGEGG=EEDFGFEFEFEEGEFFFFEFAD:BEDGBGFDF,FFAF7<FF;(GEGFFD8?FBF?/1-GGD8'FF'ECGA:AF<FF&AD\n"
              + "@CL100031435L2C001R001_63565 :TGGTCCTTGC\n"
              + "AAAGCTCTAGGGGAGGCCGACCTCTCAGCTTTTGGAGTCGGTGATACAGAGGA\n"
              + "+\n"
              + "EFFEDBCFFDAC7C@DBEF6DDFDFE;EEFFF:DA=DFDE8EC;DF@E0=AEF\n"
              + "@CL100031435L2C001R001_63596 :TGGTCCTTGC\n"
              + "CCTTCTTTTTGATCCTACAAGTGAAGCCATGTGGAAATGATTGAGAGTGACGTCACGAGTTCAGTTGTGGGGGGCGGCGGCCGGGTTGTCCGGTCGGGAC\n"
              + "+\n"
              + "CBEC>=EFCED@?CDDF<EFA>4CDD8E8AAC>:4FDD>=CDD4F7CEF;CAB=FB>:CAA?F,EFE?;&@D7>5:?/==(C>.<&26'-=?8%-+E42F",
      	output_dir.toPath(), true, "fastq-10", "gz" );
        
        FastqScanner fs = new MyScanner( expected_reads );
        RawReadsFile rf1 = new RawReadsFile();
        rf1.setFilename( f1.toFile().getCanonicalPath() );

        ValidationResult vr = new ValidationResult();

        fs.checkFiles( vr, rf1 );

        Assert.assertFalse( fs.getPaired() );
        Assert.assertEquals( 0, vr.count(Severity.ERROR) );
    }

    
    @Test public void 
    testPairWithDuplication() throws Throwable
    {
        URL  url2 = FastqScannerTest.class.getClassLoader().getResource( "uk/ac/ebi/ena/webin/cli/rawreads/EP0_GTTCCTT_S1.txt.dup.gz" );
        URL  url1 = FastqScannerTest.class.getClassLoader().getResource( "uk/ac/ebi/ena/webin/cli/rawreads/EP0_GTTCCTT_S2.txt.gz" );

        FastqScanner fs = new MyScanner( expected_reads );
        
        RawReadsFile rf1 = new RawReadsFile();
        rf1.setFilename( new File( url1.getFile() ).getCanonicalPath() );
        
        RawReadsFile rf2 = new RawReadsFile();
        rf2.setFilename( new File( url2.getFile() ).getCanonicalPath() );

        ValidationResult vr = new ValidationResult();

        fs.checkFiles( vr, rf1, rf2 );

        Assert.assertEquals( 2, vr.count(Severity.ERROR) );
    }
    

    @Test public void 
    testPairWithDuplication2() throws Throwable
    {
        URL  url2 = FastqScannerTest.class.getClassLoader().getResource( "uk/ac/ebi/ena/webin/cli/rawreads/EP0_GTTCCTT_S1.txt.dup.gz" );
        //URL  url1 = GenomeAssemblyWebinCliTest.class.getClassLoader().getResource( "uk/ac/ebi/ena/webin/cli/rawreads/EP0_GTTCCTT_S2.txt.gz" );
        FastqScanner fs = new MyScanner( expected_reads );
        //RawReadsFile rf1 = new RawReadsFile();
        //rf1.setFilename( new File( url1.getFile() ).getCanonicalPath() );
        
        RawReadsFile rf2 = new RawReadsFile();
        rf2.setFilename( new File( url2.getFile() ).getCanonicalPath() );

        ValidationResult vr = new ValidationResult();

        fs.checkFiles( vr, rf2, rf2 );

        Assert.assertEquals( 2, vr.count(Severity.ERROR) );
    }
    
    
    private Path
    generateRandomFastq(int number_of_reads,
                        int min_name_len,
                        int max_name_len,
                        int read_len) throws IOException
    {
        Path result = Files.createTempFile( "TEMP", ".fastq" );
        
        while( number_of_reads--> 0 )
        {
            StringBuilder read = new StringBuilder();

            read.append( "@" )
                .append( ThreadLocalRandom.current()
                                          .ints( ThreadLocalRandom.current().nextInt( min_name_len, max_name_len ), 55, 95 )
                                          .mapToObj( e -> Character.toString((char) e))
                                          .collect( Collectors.joining() ) )
                .append( '\n' )
                .append( ThreadLocalRandom.current()
                        .ints( read_len, 0, 5 )
                        .mapToObj( e -> e == 0 ? "A" : e == 1 ? "C" : e == 2 ? "G" : e == 3 ? "T" : "N" )
                        .collect( Collectors.joining() ) )
                .append( '\n' )
                .append( '+' )
                .append( '\n' )
                .append( ThreadLocalRandom.current()
                        .ints( read_len, 33, 33 + 64 )
                        .mapToObj( e -> Character.toString((char) e))
                        .collect( Collectors.joining() ) )
                .append( '\n' );
            Files.write( result, read.toString().getBytes(), StandardOpenOption.SYNC, StandardOpenOption.APPEND );
        }        
        return result;
            
    }

    //TODO remove probabilistic nature
    @Test public void 
    testGeneratedSingleDuplications() throws Throwable
    {
        FastqScanner fs = new MyScanner( expected_reads );
        RawReadsFile rf = new RawReadsFile();
        Path path = generateRandomFastq( 1000, 2, 3, 80 );
        rf.setFilename( path.toString() );

        ValidationResult vr = new ValidationResult();

        fs.checkFiles( vr, rf );

        Assert.assertNotEquals( 1, vr.count(Severity.ERROR) );
    }

    @Test
    public void testInvalid() throws Throwable {
        URL  url1 = FastqScannerTest.class.getClassLoader().getResource( "invalid.fastq.gz" );

        RawReadsFile rf = new RawReadsFile();

        rf.setFilename( new File( url1.getFile() ).getCanonicalPath() );

        ValidationResult vr = new ValidationResult();

        FastqScanner fastqScanner = new MyScanner( expected_reads );
        fastqScanner.checkFiles( vr, rf );

        Assert.assertFalse(vr.isValid());
    }

    @Test
    public void testSinglePairedWithUnpairedFastq() throws Throwable {
        File output_dir = createOutputFolder();

        Path f1 = saveRandomized(
                "@paired/1\n" +
                        "ACGT\n" +
                        "+\n" +
                        "1234\n" +
                        "@paired/2\n" +
                        "ACGT\n" +
                        "+\n" +
                        "1234\n" +
                        "@unpaired\n" +
                        "ACGT\n" +
                        "+\n" +
                        "1234", output_dir.toPath(), true, "fastq-1", "gz" );

        RawReadsFile rf1 = new RawReadsFile();
        rf1.setFilename( f1.toFile().getCanonicalPath() );

        ValidationResult vr = new ValidationResult();

        FastqScanner fs = new MyScanner( expected_reads );
        fs.checkFiles( vr, rf1 );

        Assert.assertFalse( fs.getPaired() );
        Assert.assertTrue(vr.isValid());
    }

    @Test
    public void test4PairedFastqs() throws Throwable {
        URL url1 = FastqScannerTest.class.getClassLoader().getResource( "10x/4fastq/I1.fastq" );
        URL url2 = FastqScannerTest.class.getClassLoader().getResource( "10x/4fastq/R1.fastq" );
        URL url3 = FastqScannerTest.class.getClassLoader().getResource( "10x/4fastq/R2.fastq" );
        URL url4 = FastqScannerTest.class.getClassLoader().getResource( "10x/4fastq/R3.fastq" );

        FastqScanner fs = new MyScanner( expected_reads );

        RawReadsFile rf1 = new RawReadsFile();
        rf1.setFilename( new File( url1.getFile() ).getCanonicalPath() );

        RawReadsFile rf2 = new RawReadsFile();
        rf2.setFilename( new File( url2.getFile() ).getCanonicalPath() );

        RawReadsFile rf3 = new RawReadsFile();
        rf3.setFilename( new File( url3.getFile() ).getCanonicalPath() );

        RawReadsFile rf4 = new RawReadsFile();
        rf4.setFilename( new File( url4.getFile() ).getCanonicalPath() );

        ValidationResult vr = new ValidationResult();

        fs.checkFiles( vr, rf1, rf2, rf3, rf4 );

        Assert.assertTrue( vr.isValid() );
    }
}
