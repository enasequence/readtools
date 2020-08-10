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
package uk.ac.ebi.ena.readtools.loader.cg.reads;

import java.io.IOException;
import java.io.InputStream;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

public class 
ReadHeader extends CompleteGenomicsBase
{
    public interface 
    PropertyName
    {
        public String ASSEMBLY_ID = "ASSEMBLY_ID"; // Name of the assembly "<assembly-name>-ASM". For example, "GS000000474-ASM".
        public String BATCH_FILE_NUMBER = "BATCH_FILE_NUMBER"; // Number of the batch of a split data file. Positive 1-based integer.
        public String BATCH_OFFSET = "BATCH_OFFSET"; // Offset of the first record in a batch to the position of the record in a non-split file. Positive 0-based integer.
        public String FIELD_SIZE = "FIELD_SIZE"; // Size of the lane fields. Positive integer.
        public String FORMAT_VERSION = "FORMAT_VERSION"; // Version number of the file format. Two or more digits separated by periods. For example, �0.6�.
        public String GENERATED_AT = "GENERATED_AT"; // Date and time of the assembly. Year-Month-Day Time. For example �2010-Sep-08 20:27:52.457773�.
        public String GENERATED_BY = "GENERATED_BY"; // Assembly pipeline component that generated the output. Alpha-numeric string.
        public String LANE = "LANE"; // Identifier of the slide lane from which the reads were extracted
        public String LIBRARY = "LIBRARY"; // Identifier of the library from which the DNBs were generated
        public String SAMPLE = "SAMPLE"; // Complete Genomics identifier of the sample from which the library was created �GSXXXXX-DNA_YZZ� 
        public String SLIDE = "SLIDE"; // Flow slide identification code
        public String SOFTWARE_VERSION = "SOFTWARE_VERSION "; // Assembly pipeline build number. Two or more digits separated by periods.
        public String TYPE = "TYPE"; // Indicates the type of data contained in the file �READS�: reads file.
    }
    
    InputStream istream;
    final char   PROPERTY_FLAG = '#';
    final char   COLUMN_ROW_FLAG = '>';
    final String SEPARATOR = "\t";
    final int    FLAG_POS = 0;
    final int    VALUES_POS = 1;
    final Map<String, String> properties = new HashMap<String, String>(); 
    String[]     columns;
    
    
    public
    ReadHeader( InputStream istream )
    {
        super( null );
        this.istream = istream;
        
    }
    
    public Object
    getKey()
    {
        return null;
    }
    
    
    public void 
    read() throws IOException
    {
cycle:  do
        {
            String line = readLine( istream ).trim();
            
            if( line.length() > 0 )
            {
                char flag = line.charAt( FLAG_POS );
                String value = line.substring( VALUES_POS );
                switch( flag )
                {
                case PROPERTY_FLAG:
                    String[] values = value.split( SEPARATOR );
                    properties.put( values[ 0 ], values[ 1 ] );
                    break;
                    
                case COLUMN_ROW_FLAG:
                    columns = value.split( SEPARATOR );
                    break cycle;
                    
                default:
                    throw new IOException();
                }
            }
        }while( true );
    }
    
    
    public String
    getFileProrerty( String key )
    {
        return properties.get( key );
    }
    
    
    public Map<String, String>
    getFileProrerties()
    {
        return properties;
    }
    
    
    public String[]
    getColumnHeaders()
    {
        return columns;
    }
    
    
    
    @Override
    public String 
    toString()
    {
        return String.format( "Properties: %s\nColumns: %s", properties, Arrays.asList( columns ) );
    }
}
