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
package uk.ac.ebi.ena.readtools.sam;

import htsjdk.samtools.util.Log;
import org.junit.Assert;
import org.junit.Test;

import java.io.File;
import java.io.IOException;
import java.net.URISyntaxException;
import java.net.URL;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;

/**
 * The expected TotalReadCount and TotalBaseCount in the following tests were calculated using these tools on sra-login:<br/>
 * - ~/tools/putils/bam_stats<br/>
 * - ~/tools/putils/cram_stats
 */
public class
Sam2FastqTest
{
	@Test
	public void convertCram2Fastq() throws Exception {
		GeneratedFastqResult genRes = generateFastqFiles("SRR2989699.cram");

		File output = genRes.output;

		Assert.assertEquals( "@SRR2989699.5\n"
						   + "CGCCACGAGCTGGTTGTCTATGGGACAAGTGATGTGGTTGATAACCTCCCATTGCTATCTCA\n"
						   + "+\n"
						   + "EEEEEBEEEEEEEADDEEDEFFBABDEDEBBDADAFEBFAFAFFBFEFFEFDFEF=FAFFFD\n"
						   + "@SRR2989699.4\n"
						   + "CGAAGAACAGATGGCTCGTACCAAGCAAACCGCTCGTAAGTCCACCGGAGGTAAAGCTCCAA\n"
						   + "+\n"
						   + "FGGFGFFGBGFEFGGFGA=GFEGGGFGGGFGGFEGEDGDFFFFFFDBDDFE?EBE?B@C@AE\n"
						   + "@SRR2989699.1\n"
						   + "TCGTTCTTTGCCTCCGGGATTGTGAGCCAACATAATGCAATCGACATGGACTACTTGTTTGG\n"
						   + "+\n"
						   + "5DEGEGAGDFGGGGFDGG?GFGFDGGGGAGGEGEGGEGGGGFGDFDFGF=F=F=FFF5FADE\n"
						   + "@SRR2989699.6\n"
						   + "CAACAACCTTCTCAGCAACTAAGAAAGCAGAGTAAAACCCAACACCAAACTGTCCGACAATG\n"
						   + "+\n"
						   + "GGGGEGFGDFFFFEFFFGEEFFEFFFGEGGGGGGFCFEFFFDGGEDFEFDEDEFEEAD=ECD\n",
				             new String( Files.readAllBytes( Paths.get( output.getPath() + ".fastq" ) ), StandardCharsets.UTF_8 ) );

		Assert.assertEquals(4l, genRes.sam2Fastq.getTotalReadCount());
		Assert.assertEquals(248l, genRes.sam2Fastq.getTotalBaseCount());
	}

	@Test
	public void convertCram2Fastq2Files() throws Exception {
		String baseDir = "cram2fastq/2fastq/";
		String fileNamePrefix = "28239_1822";
		String fileExt = ".cram";

		GeneratedFastqResult genRes = generateFastqFiles(baseDir + fileNamePrefix + fileExt);

		File output = genRes.output;

		assertFastqResult(baseDir + fileNamePrefix, output.getPath(), "_1");
		assertFastqResult(baseDir + fileNamePrefix, output.getPath(), "_2");

		Assert.assertEquals(6l, genRes.sam2Fastq.getTotalReadCount());
		Assert.assertEquals(900l, genRes.sam2Fastq.getTotalBaseCount());
	}

	@Test
	public void convertBam2Fastq1File() throws Exception {
		String baseDir = "bam2fastq/1fastq/";
		String fileNamePrefix = "S0567a_E1_L1__aln.sort.mapped.rmdupse_adna_v2";
		String fileExt = ".bam";

		GeneratedFastqResult genRes = generateFastqFiles(baseDir + fileNamePrefix + fileExt);

		File output = genRes.output;

		assertFastqResult(baseDir + fileNamePrefix, output.getPath(), "");

		Assert.assertEquals(20l, genRes.sam2Fastq.getTotalReadCount());
		Assert.assertEquals(1249l, genRes.sam2Fastq.getTotalBaseCount());
	}

	@Test
	public void convertBam2Fastq2Files() throws Exception {
		String baseDir = "bam2fastq/2fastq/";
		String fileNamePrefix = "8855_124";
		String fileExt = ".bam";

		GeneratedFastqResult genRes = generateFastqFiles(baseDir + fileNamePrefix + fileExt);

		File output = genRes.output;

		assertFastqResult(baseDir + fileNamePrefix, output.getPath(), "_1");
		assertFastqResult(baseDir + fileNamePrefix, output.getPath(), "_2");

		Assert.assertEquals(234l, genRes.sam2Fastq.getTotalReadCount());
		Assert.assertEquals(7020l, genRes.sam2Fastq.getTotalBaseCount());
	}

	@Test
	public void convertBam2Fastq3Files() throws Exception {
		String baseDir = "bam2fastq/3fastq/";
		String fileNamePrefix = "M2241_BLV_sense";
		String fileExt = ".bam";

		GeneratedFastqResult genRes = generateFastqFiles(baseDir + fileNamePrefix + fileExt);

		File output = genRes.output;

		assertFastqResult(baseDir + fileNamePrefix, output.getPath(), "_1");
		assertFastqResult(baseDir + fileNamePrefix, output.getPath(), "_2");
		assertFastqResult(baseDir + fileNamePrefix, output.getPath(), "");

		Assert.assertEquals(129l, genRes.sam2Fastq.getTotalReadCount());
		Assert.assertEquals(13029l, genRes.sam2Fastq.getTotalBaseCount());
	}

	@Test
	public void convertSam2FastqReversePairedReads() throws Exception {
		String baseDir = "sam2fastq/";
		String fileNamePrefix = "reverse-read";
		String fileExt = ".sam";

		File output = generateFastqFiles(baseDir + fileNamePrefix + fileExt).output;

		assertFastqResult(baseDir + fileNamePrefix, output.getPath(), "_1");
		assertFastqResult(baseDir + fileNamePrefix, output.getPath(), "_2");
	}

	private GeneratedFastqResult generateFastqFiles(String source) throws Exception {
		File output = File.createTempFile( "FASTQ", "FASTQ" );
		output.delete();

		URL url = Sam2FastqTest.class.getClassLoader().getResource( source);
		File file = new File( url.getFile() );

		Sam2Fastq.Params params = new Sam2Fastq.Params();
		params.samFile = file;
		params.reverse = true;
		params.nofStreams = 3;
		params.fastqBaseName = output.getPath();

		Log.setGlobalLogLevel( params.logLevel );

		Sam2Fastq sam2Fastq = new Sam2Fastq();
		sam2Fastq.create( params );

		return new GeneratedFastqResult(sam2Fastq, output);
	}

	private void assertFastqResult(String expectedFilePrefix, String actualFilePrefix, String fileIndex) throws URISyntaxException, IOException {
		String expectedFastq = new String(Files.readAllBytes(Paths.get(Sam2FastqTest.class.getClassLoader()
				.getResource( expectedFilePrefix + fileIndex + ".fastq" ).toURI())), StandardCharsets.UTF_8);

		String actualFastq = new String(Files.readAllBytes(Paths.get( actualFilePrefix + fileIndex + ".fastq" ) ),StandardCharsets.UTF_8);

		Assert.assertEquals(expectedFastq, actualFastq);
	}

	private static class GeneratedFastqResult {
		public Sam2Fastq sam2Fastq;
		public File output;

		public GeneratedFastqResult(Sam2Fastq sam2Fastq, File output) {
			this.sam2Fastq = sam2Fastq;
			this.output = output;
		}
	}
}
