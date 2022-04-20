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

import uk.ac.ebi.ena.readtools.common.reads.normalizers.htsjdk.StandardQualityNormalizer;
import uk.ac.ebi.ena.readtools.loader.common.InvalidBaseCharacterException;

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
}
