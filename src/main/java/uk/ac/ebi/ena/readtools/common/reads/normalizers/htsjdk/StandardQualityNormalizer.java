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
package uk.ac.ebi.ena.readtools.common.reads.normalizers.htsjdk;

import htsjdk.samtools.SAMUtils;

import uk.ac.ebi.ena.readtools.common.reads.QualityNormalizer;

/**
 * Converts printable qualities in Sanger fastq format to binary phred scores.
 */
public class StandardQualityNormalizer implements QualityNormalizer {

    @Override
    public void normalize(byte[] qualities) {
        SAMUtils.fastqToPhred(qualities);
    }
}
