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

public interface Converter {
    int YIELD_CYCLES_FOR_ERROR_CHECKING = 362; // 16384;

    /**
     * Get the total number of reads that were read.
     *
     * @return
     */
    long getReadCount();

    /**
     * Get the total number of bases that were read.
     *
     * @return
     */
    long getBaseCount();

    void run();
    void runOnce();
    boolean isDone();
}
