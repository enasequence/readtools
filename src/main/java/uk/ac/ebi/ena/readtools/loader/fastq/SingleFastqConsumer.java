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
SingleFastqConsumer implements DataConsumer<DataSpot, IlluminaSpot>
{
    DataConsumer<IlluminaSpot, ?> dataConsumer;
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
        
        IlluminaSpot ispot = IlluminaSpot.initSingle();
        int slash_idx = spot.bname.lastIndexOf( '/' );
        ispot.name = slash_idx == -1 ? spot.bname 
                                     : spot.bname.substring( 0, slash_idx );
        ispot.bases = spot.bases;
        ispot.quals = spot.quals;
        ispot.read_start[ IlluminaSpot.FORWARD ] = 0;
        ispot.read_length[ IlluminaSpot.FORWARD ] = spot.bases.length();
        ispot.read_name[ IlluminaSpot.FORWARD ] = spot.bname;
        if( null != dataConsumer)
            dataConsumer.consume( ispot );
        else
            System.out.println( ispot );
    }


    @Override
    public void
    setConsumer(DataConsumer<IlluminaSpot, ?> dataConsumer)
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
