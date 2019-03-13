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
