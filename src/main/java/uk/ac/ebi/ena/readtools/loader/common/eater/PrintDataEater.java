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
package uk.ac.ebi.ena.readtools.loader.common.eater;



public class 
PrintDataEater<T1, T2> implements DataConsumer<T1, T2>
{
    volatile boolean lock;
    
    @Override
    public void 
    cascadeErrors() throws DataConsumerException
    {
        ;
    }

 
    public boolean
    setLock()
    {
        if( lock )
            return false;
        lock = true;
        return lock;
    }
    

    public void
    resetLock()
    {
        lock = false;
    }

    
    @Override public void
    consume(T1 object ) throws DataConsumerException
    {
        System.out.println( object ); 
    }

    
    @Override public void
    setConsumer(DataConsumer<T2, ?> dataConsumer)
    {
        throw new UnsupportedOperationException();
    }


    @Override public boolean 
    isOk()
    {
        return true;
    }

}
