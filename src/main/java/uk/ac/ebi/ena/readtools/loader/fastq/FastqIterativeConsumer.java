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

import java.io.File;
import java.util.Iterator;

import uk.ac.ebi.ena.readtools.loader.common.QualityNormalizer;

public class
FastqIterativeConsumer implements Iterable<FastqSpot>
{
    public enum
    READ_TYPE
    {
      SINGLE,
      PAIRED
    };
    
    
    private File[] files;
    private int spill_page_size = 5000000;
    private File tmp_folder = new File( "." );
    private READ_TYPE read_type;
    private QualityNormalizer[] normalizers;
    
    
    @Override
    public Iterator<FastqSpot>
    iterator()
    {
        try
        {
            return (Iterator<FastqSpot>) new FastqIterativeConsumerIterator( tmp_folder,
                                                                        spill_page_size, 
                                                                        read_type, 
                                                                        files,
                                                                        normalizers );
        } catch( Throwable t )
        {
            t.printStackTrace();
        }
        
        return null;
    }


    public void
    setFiles( File files[] )
    {
        this.files = files;
    }

    
    public File[]
    getFiles()
    {
        return this.files;
    }
    
    
    public int 
    getSpillPageSize()
    {
        return spill_page_size;
    }


    public void 
    setSpillPageSize( int spill_page_size )
    {
        this.spill_page_size = spill_page_size;
    }


    public File 
    getTmpFolder()
    {
        return tmp_folder;
    }


    public void 
    setTmpFolder( File tmp_folder )
    {
        this.tmp_folder = tmp_folder;
    }


    public READ_TYPE 
    getReadType()
    {
        return read_type;
    }


    public void 
    setReadType( READ_TYPE read_type )
    {
        this.read_type = read_type;
    }


    public QualityNormalizer[] 
    getNormalizers()
    {
        return normalizers;
    }


    public void 
    setNormalizers( QualityNormalizer[] normalizers )
    {
        this.normalizers = normalizers;
    }
}
