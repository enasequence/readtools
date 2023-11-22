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
package uk.ac.ebi.ena.readtools.loader.common.converter;

import java.io.ByteArrayInputStream;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;

import org.junit.Assert;
import org.junit.Test;

import uk.ac.ebi.ena.readtools.common.reads.QualityNormalizer;
import uk.ac.ebi.ena.readtools.common.reads.normalizers.htsjdk.StandardQualityNormalizer;
import uk.ac.ebi.ena.readtools.loader.common.InvalidBaseCharacterException;
import uk.ac.ebi.ena.readtools.loader.fastq.Read;

public class
ReadReaderTest {

	@Test
	public void testCASAVA1_8() {
		//CASAVA 1.8
		Assert.assertTrue( ReadReader.p_casava_1_8_name.matcher( "@EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG" ).matches() );
		Assert.assertTrue( ReadReader.p_casava_1_8_name.matcher( "@EAS139:136:FC706VJ:2:2104:15343:197393 2:Y:18:ATCACG" ).matches() );

		//extended version of CASAVA 1.8
		Assert.assertTrue( ReadReader.p_casava_1_8_name.matcher( "@M00825:71:000000000-AARLA:1:1101:16089:1603 1:N:0:331" ).matches() );
		Assert.assertTrue( ReadReader.p_casava_1_8_name.matcher( "@M00825:71:000000000-AARLA:1:1101:16089:1603 2:N:0:331" ).matches() );

		Assert.assertTrue( ReadReader.p_casava_1_8_name.matcher( "@M00825:71:000000000-AARLA:1:1101:16089:1603     2:N:0:331" ).matches() );
		Assert.assertTrue( ReadReader.p_casava_1_8_name.matcher( "@M00825:71:000000000-AARLA:1:1101:16089:1603\t2:N:0:331" ).matches() );
		Assert.assertTrue( ReadReader.p_casava_1_8_name.matcher( "@M00825:71:000000000-AARLA:1:1101:16089:1603\t\t2:N:0:331" ).matches() );

		Assert.assertTrue( ReadReader.p_casava_1_8_name.matcher( "@M00825:71:000000000-AARLA:1:1101:16089:1603 1:N:0:331 COMMENT" ).matches() );

		Assert.assertTrue( ReadReader.p_base_name
				.matcher( "@A00372:119:HNJM2DMXX:1:1101:17282:5055:N:0:CGGTTACGGC+AAGACTATAG#0/1" ).matches() );

		Assert.assertTrue( ReadReader.p_casava_1_8_name
				.matcher( "@EAS139:136:FC706VJ:2:5:1000:12850 1:Y:18:ATCACGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAABBBBCCCC?<A?BC?7@@???????DBBA@@@@A@@" ).matches() );


		Assert.assertFalse( ReadReader.p_casava_1_8_name
				.matcher( "@A00372:119:HNJM2DMXX:1:1101:9588:1752:N:0:CGGTTACGGC+AAGACTATAG#0/1" ).matches() );
		Assert.assertFalse( ReadReader.p_casava_1_8_name
				.matcher( "@A00372:119:HNJM2DMXX:1:1101:17282:5055:N:0:CGGTTACGGC+AAGACTATAG#0/1" ).matches() );

		Assert.assertFalse( ReadReader.p_casava_1_8_name
				.matcher( "@ 9:N:0:ACAGCAAC" ).matches() );

		Assert.assertTrue( ReadReader.p_casava_1_8_name
				.matcher( "@A00730:546:HWCTCDRXY:2:2101:1090:1031 9:N:0:ACAGCAAC" ).matches() );
		Assert.assertTrue( ReadReader.p_casava_1_8_name
				.matcher( "@E00528:414:HVNJLCCXY:1:1101:7598:1854 1:N:0" ).matches() );
		Assert.assertTrue( ReadReader.p_casava_1_8_name
				.matcher( "@A00730:546:HWCTCDRXY:2:2101:1090:1031:ACAGCAAC+ACAGCAAC 4:N:0:ACAGCAAC" ).matches() );
		Assert.assertTrue( ReadReader.p_casava_1_8_name
				.matcher( "@A00730:546:HWCTCDRXY:2:2101:1090:1031:ACAGCAAC+ACAGCAAC 422:N:0:ACAGCAAC" ).matches() );
		Assert.assertTrue( ReadReader.p_casava_1_8_name
				.matcher( "@A00730:546:HWCTCDRXY:2:2101:1090:1031:ACAGCAAC+ACAGCAAC\t422:N:0:ACAGCAAC" ).matches() );
		Assert.assertFalse( ReadReader.p_casava_1_8_name
				.matcher( "@A00730:546:HWCTCDRXY:2:2101:1090:1031:ACAGCAAC+ACAGCAAC\n422:N:0:ACAGCAAC" ).matches() );

		Assert.assertTrue( ReadReader.p_casava_1_8_name
				.matcher( "@JAXVAFT.MTP3.D21.48_239647 M00990:616:000000000-JKYVV:1:1101:20204:2129 1:N:0:TTACTGTGCG+GATTCCTA" ).matches() );

		Assert.assertFalse( ReadReader.p_casava_1_8_name
				.matcher( "@A00953:544:HMTFHDSX3:2:1101:23981:1814" ).matches() );
		Assert.assertTrue( ReadReader.p_base_name
				.matcher( "@A00953:544:HMTFHDSX3:2:1101:23981:1814" ).matches() );
	}

	@Test
	public void testQualityScoreNormalization() throws Exception {
		QualityNormalizer normalizer = new StandardQualityNormalizer();
		String qualityScores = "FFFFFFFFFFFFFFFF";
		String expectedNormalizedQualityScores = "%%%%%%%%%%%%%%%%";

		String fullInput = String.format("%s\n%s\n%s\n%s",
			"@RN-001",
			"AGCTAGCTAGCTAGCT",
			"+",
			qualityScores);

		try(InputStream is = new ByteArrayInputStream((fullInput).getBytes(StandardCharsets.UTF_8))) {
			ReadReader ds = new ReadReader(normalizer);

			Read read = ds.read(is);

			Assert.assertEquals(expectedNormalizedQualityScores, read.getQualityScores());
		}
	}

	@Test
	public void testInvalidBaseCharactersException() throws Exception {
		String basesInput = "AGCTUXagctuAGCTUagctux";

		String fullInput = String.format("%s\n%s\n%s\n%s",
				"@RN-001",
				basesInput,
				"+",
				"FFFFFFFFFFFFFFFFFFFF");

		InvalidBaseCharacterException ex = null;

		//DataSpot keeps reading until it sees '+' symbol on a new line which is in compliance with Fastq standard.
		try(InputStream is = new ByteArrayInputStream((fullInput).getBytes(StandardCharsets.UTF_8))) {
			ReadReader ds = new ReadReader(new StandardQualityNormalizer());
			ds.read(is);
		} catch (InvalidBaseCharacterException e) {
			ex = e;
		} catch (Exception e) {
			throw e;
		}

		Assert.assertEquals(basesInput, ex.getBases());
		Assert.assertArrayEquals(new Character[]{'X', 'x'}, ex.getInvalidCharacters().toArray());
	}

	@Test
	public void testReadNameLengthExceeds256() throws Exception {
		String basesInput = "AGCTUagctuAGCTUagctu";
		String fullInput = String.format("%s\n%s\n%s\n%s",
				"@0000000000000000000000000000000000000000000000000000000000000000" +
						"0000000000000000000000000000000000000000000000000000000000000000" +
						"0000000000000000000000000000000000000000000000000000000000000000" +
						"0000000000000000000000000000000000000000000000000000000000000000" +
						"1",
				basesInput,
				"+",
				"FFFFFFFFFFFFFFFFFFFF");

		//DataSpot keeps reading until it sees '+' symbol on a new line which is in compliance with Fastq standard.
		try(InputStream is = new ByteArrayInputStream((fullInput).getBytes(StandardCharsets.UTF_8))) {
			ReadReader ds = new ReadReader(new StandardQualityNormalizer());
			ds.read(is);
			Assert.fail();
		} catch (ConverterException e) {
			Assert.assertTrue(e.getMessage().contains("Line's length exceeds 256 characters"));
		}
	}

	@Test
	public void testReadNameLength256() throws Exception {
		String basesInput = "AGCTUagctuAGCTUagctu";
		String fullInput = String.format("%s\n%s\n%s\n%s",
				"@0000000000000000000000000000000000000000000000000000000000000000" +
						"0000000000000000000000000000000000000000000000000000000000000000" +
						"0000000000000000000000000000000000000000000000000000000000000000" +
						"0000000000000000000000000000000000000000000000000000000000000000",
				basesInput,
				"+",
				"FFFFFFFFFFFFFFFFFFFF");

		//DataSpot keeps reading until it sees '+' symbol on a new line which is in compliance with Fastq standard.
		try(InputStream is = new ByteArrayInputStream((fullInput).getBytes(StandardCharsets.UTF_8))) {
			ReadReader ds = new ReadReader(new StandardQualityNormalizer());
			ds.read(is);
		}
	}
}
