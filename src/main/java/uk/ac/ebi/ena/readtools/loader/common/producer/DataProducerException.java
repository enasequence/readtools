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
package uk.ac.ebi.ena.readtools.loader.common.producer;

public class
DataProducerException extends RuntimeException
{
    private static final long serialVersionUID = 1L;
    private long   line_no = -1;
    
    public DataProducerException(long line_no, String value )
    {
        super( value );
        this.line_no = line_no;
    }

    public DataProducerException(long line_no )
    {
        super();
        this.line_no = line_no;
    }

    
    public DataProducerException(long line_no, Throwable cause )
    {
        super( cause );
        this.line_no = line_no;
    }


    public DataProducerException(String message, Throwable cause )
    {
        super( message, cause );
    }
    
    
    public DataProducerException(Throwable cause )
    {
        super( cause );
    }
   
    
    public String
    toString()
    {
        return line_no < 0 ? super.toString() 
                           : String.format( "Line No. : %d - %s", line_no, super.toString() );
    }
    
    
    public long
    getLineNo()
    {
    	return this.line_no;
    }
}
