/*******************************************************************************
 * Copyright 2013 EMBL-EBI
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *   http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 ******************************************************************************/
package htsjdk.samtools.cram.encoding.reader;

public abstract class AbstractFastqReader{
	public boolean reverseNegativeReads = true;
	public boolean appendSegmentIndexToReadNames = true;

	public int readLength;
	public byte[] readName;

	public static final int maxReadBufferLength = 1024 * 1024;
	public byte[] bases = new byte[maxReadBufferLength];
	public byte[] scores = new byte[maxReadBufferLength];

	public int defaultQS = '?';

	/**
	 * For now this is to identify the right buffer to use.
	 * 
	 * @param flags
	 *            read bit flags
	 * @return 0 for non-paired or other rubbish which could not be reliably
	 *         paired, 1 for first in pair and 2 for second in pair
	 */
	protected int getSegmentIndexInTemplate(int flags) {
		if ((flags & 1) == 0)
			return 0;

		if ((flags & 64) != 0)
			return 1;
		else
			return 2;
	}

	/**
	 * Write the read. The read here is basically a fastq read with an addition
	 * of SAM bit flags. Specific implementations should take care of further
	 * cashing/pairing/filtering and actual writing of reads.
	 * <p>
	 * The contract is
	 * <ul>
	 * <li>no supplementary reads will appear in this method.
	 * <li>reads on negative strand will be reverse complimented to appear as if
	 * on positive strand.
	 * </ul>
	 * 
	 * @param name
	 *            read name
	 * @param flags
	 *            SAM bit flags
	 * @param bases
	 *            read bases
	 * @param scores
	 *            fastq quality scores (phred+33)
	 */
	public abstract void writeRead(byte[] name, int flags, byte[] bases, byte[] scores);

	public abstract void finish();
}
