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
package uk.ac.ebi.ena.readtools.fastq.enafastq2sam;

import org.apache.commons.lang.NotImplementedException;
import uk.ac.ebi.ena.readtools.loader.common.QualityConverter;
import uk.ac.ebi.ena.readtools.loader.common.QualityNormalizer;

public enum ILLUMINA_SETTINGS {
    SANGER,
    SOLEXA,
    ILLUMINA_1_3,
    ILLUMINA_1_5;
    
    public static ILLUMINA_SETTINGS valueOf( QualityNormalizer normalizer, QualityConverter converter ) {
        switch( normalizer ) {
            case SANGER:
                return SANGER;

            case SOLEXA:
                return SOLEXA;

            case ILLUMINA_1_3:
                return ILLUMINA_1_3;

            case ILLUMINA_1_5:
                return ILLUMINA_1_5;

            case X:
                switch( converter ) {
                    case PHRED:
                        return SANGER;

                    default:
                        throw new NotImplementedException();
                }

            case X_2:
                switch( converter ) {
                    case PHRED:
                        return ILLUMINA_1_3;

                    case LOGODDS:
                        return SOLEXA;

                    default:
                        throw new NotImplementedException();
                }

            default:
                throw new NotImplementedException();
        }
    }
    
    public QualityNormalizer getNormalizer() {
        switch( this )
        {
            case SANGER:
                return QualityNormalizer.SANGER;
                
            case SOLEXA:
                return QualityNormalizer.SOLEXA;
            
            case ILLUMINA_1_3:
                return QualityNormalizer.ILLUMINA_1_3;
            
            case ILLUMINA_1_5:
                return QualityNormalizer.ILLUMINA_1_5;
               
            default:
                throw new NotImplementedException();
        }
    }
}
