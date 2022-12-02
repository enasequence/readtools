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
package uk.ac.ebi.ena.readtools.loader.fastq;


import uk.ac.ebi.ena.readtools.loader.common.writer.ReadWriter;
import uk.ac.ebi.ena.readtools.loader.common.writer.ReadWriterException;
import uk.ac.ebi.ena.readtools.loader.common.writer.Spot;

public class
PrintReadWriter<T1 extends Spot, T2 extends Spot> implements ReadWriter<T1, T2>
{
    @Override
    public void 
    cascadeErrors() throws ReadWriterException
    {
    }

    @Override public void
    write(T1 spot ) throws ReadWriterException
    {
        System.out.println( spot );
    }

    
    @Override public void
    setWriter(ReadWriter<T2, ?> readWriter)
    {
        throw new UnsupportedOperationException();
    }
}
