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

import uk.ac.ebi.ena.readtools.loader.common.consumer.DataConsumer;
import uk.ac.ebi.ena.readtools.loader.common.consumer.DataConsumerException;

public class
SingleFastqConsumer implements DataConsumer<DataSpot, FastqSpot>
{
    DataConsumer<FastqSpot, ?> dataConsumer;
    boolean is_ok = true;
    
    @Override
    public void 
    cascadeErrors() throws DataConsumerException
    {
        //empty, no errors here
        
        if( null != dataConsumer)
            dataConsumer.cascadeErrors();
    }

    
    @Override
    public void
    consume(DataSpot spot ) throws DataConsumerException
    {
        int slash_idx = spot.name.lastIndexOf( '/' );

        FastqSpot fastqSpot = new FastqSpot(slash_idx == -1 ? spot.name : spot.name.substring( 0, slash_idx ), spot);

        if( null != dataConsumer)
            dataConsumer.consume( fastqSpot );
        else
            System.out.println( fastqSpot );
    }


    @Override
    public void
    setConsumer(DataConsumer<FastqSpot, ?> dataConsumer)
    {
        this.dataConsumer = dataConsumer;
    }


    @Override
    public boolean 
    isOk()
    {
        return null == dataConsumer ? is_ok : is_ok && dataConsumer.isOk();
    }
}
