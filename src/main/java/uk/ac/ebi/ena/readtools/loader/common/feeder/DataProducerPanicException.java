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
package uk.ac.ebi.ena.readtools.loader.common.feeder;

public class
DataProducerPanicException extends DataProducerException
{
    private static final long serialVersionUID = 1L;
    private final String thread_name;
    
    
    public DataProducerPanicException(String value )
    {
        super( -1, value );
        thread_name = Thread.currentThread().getName();
    }

    
    public DataProducerPanicException()
    {
        super( -1 );
        thread_name = Thread.currentThread().getName();
    }
    
    
    public DataProducerPanicException(Throwable cause )
    {
        super( cause );
        thread_name = Thread.currentThread().getName();
    }
    
    
    public DataProducerPanicException(String message, Throwable cause )
    {
        super( message, cause );
        thread_name = Thread.currentThread().getName();
    }
    
    
    public String
    toString()
    {
        return String.format( "%s: %s", thread_name, super.toString() ); 
    }
    
}
