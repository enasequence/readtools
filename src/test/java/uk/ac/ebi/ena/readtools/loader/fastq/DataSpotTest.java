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
package uk.ac.ebi.ena.readtools.loader.fastq;

import org.junit.Assert;
import org.junit.Test;


public class
DataSpotTest 
{
	@Test public void
	testCASAVA1_8()
	{
		//CASAVA 1.8
		Assert.assertTrue( DataSpot.p_casava_1_8_name.matcher( "@EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG" ).matches() );
		Assert.assertTrue( DataSpot.p_casava_1_8_name.matcher( "@EAS139:136:FC706VJ:2:2104:15343:197393 2:Y:18:ATCACG" ).matches() );

		//extended version of CASAVA 1.8
		Assert.assertTrue( DataSpot.p_casava_1_8_name.matcher( "@M00825:71:000000000-AARLA:1:1101:16089:1603 1:N:0:331" ).matches() );
		Assert.assertTrue( DataSpot.p_casava_1_8_name.matcher( "@M00825:71:000000000-AARLA:1:1101:16089:1603 2:N:0:331" ).matches() );
		
		Assert.assertTrue( DataSpot.p_casava_1_8_name.matcher( "@M00825:71:000000000-AARLA:1:1101:16089:1603 1:N:0:331 COMMENT" ).matches() );
	}
	
}
