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

public class
SingleFastqConsumer implements ReadWriter<Read, PairedRead>
{
    ReadWriter<PairedRead, ?> readWriter;
    boolean is_ok = true;
    
    @Override
    public void 
    cascadeErrors() throws ReadWriterException
    {
        //empty, no errors here
        
        if( null != readWriter)
            readWriter.cascadeErrors();
    }

    
    @Override
    public void
    write(Read spot ) throws ReadWriterException
    {
        PairedRead pairedRead = new PairedRead(spot.key, spot);

        if( null != readWriter)
            readWriter.write(pairedRead);
        else
            System.out.println(pairedRead);
    }


    @Override
    public void
    setWriter(ReadWriter<PairedRead, ?> readWriter)
    {
        this.readWriter = readWriter;
    }


    @Override
    public boolean 
    isOk()
    {
        return null == readWriter ? is_ok : is_ok && readWriter.isOk();
    }
}
