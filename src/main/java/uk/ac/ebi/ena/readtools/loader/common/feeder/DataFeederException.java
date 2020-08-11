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
DataFeederException extends Exception
{
    private static final long serialVersionUID = 1L;
    private long   line_no = -1;
    private final  String thread_name = Thread.currentThread().getName();
    
    public 
    DataFeederException( long line_no, String value )
    {
        super( value );
        this.line_no = line_no;
    }

    public 
    DataFeederException( long line_no )
    {
        super();
        this.line_no = line_no;
    }

    
    public 
    DataFeederException( long line_no, Throwable cause )
    {
        super( cause );
        this.line_no = line_no;
    }


    public 
    DataFeederException( String message, Throwable cause )
    {
        super( message, cause );
    }
    
    
    public 
    DataFeederException( Throwable cause )
    {
        super( cause );
    }
   
    
    public String
    toString()
    {
        return line_no < 0 ? super.toString() 
                           : String.format( "%s:%d %s", thread_name, line_no, super.toString() );
    }
    
    
    public long
    getLineNo()
    {
    	return this.line_no;
    }
}
