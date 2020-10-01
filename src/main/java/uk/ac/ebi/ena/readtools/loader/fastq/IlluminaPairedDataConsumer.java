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
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import uk.ac.ebi.ena.readtools.loader.common.eater.AbstractPagedDataConsumer;
import uk.ac.ebi.ena.readtools.loader.common.eater.DataConsumerException;


public class
IlluminaPairedDataConsumer extends AbstractPagedDataConsumer<DataSpot, IlluminaSpot>
{
    // Provided readname structure is @{readkey}{separator:1(.|/|:|_)}{index:1(0:1:2)}
        
    static final Pattern split_regexp = Pattern.compile( "^(.*)(?:[\\.|:|/|_])([12])$" ); 
    static public final int KEY = 1;
    static public final int INDEX = 2;
    
    
    public IlluminaPairedDataConsumer(File tmp_root, int spill_page_size )
    {
        super( tmp_root, spill_page_size );
    }
    
    
    public static String 
    getReadnamePart( String readname, int group ) throws DataConsumerException
    {
        Matcher m = split_regexp.matcher( readname );
        if( m.find() )
            return m.group( group );
        
        throw new DataConsumerException( String.format( "Readname [%s] does not match regexp", readname ) );
        
    }
    
    
    @Override public Object 
    getKey( DataSpot object )
    {
        try
        {   return getReadnamePart( object.bname, KEY );
        } catch (DataConsumerException de )
        {
            return object.bname;
        }
    }

    
    public List<DataSpot>
    newListBucket()
    {
        List<DataSpot> list = super.newListBucket();
        list.add( null );
        list.add( null );
        return list;
    }
    
    
    public void
    append( List<DataSpot> list, DataSpot object ) throws DataConsumerException
    {
        String index_part;
        try
        {
            index_part = getReadnamePart( object.bname, INDEX );
        } catch ( DataConsumerException de )
        {
            index_part = object.getStreamKey();
        }
        
        int index = Integer.parseInt( index_part ) - 1;
        if( null == list.get( index ) )
            list.set( index, object );
        else
            throw new RuntimeException( "Got same spot twice: " + object );
    }
    
    
    @Override
    public IlluminaSpot 
    assemble( final Object key, List<DataSpot> list ) throws DataConsumerException
    {
        IlluminaSpot i_spot = IlluminaSpot.initPaired();
        
        i_spot.read_start[ IlluminaSpot.FORWARD ] = 0;
        i_spot.read_start[ IlluminaSpot.REVERSE ] = 0;

        i_spot.read_length[ IlluminaSpot.FORWARD ] = 0;
        i_spot.read_length[ IlluminaSpot.REVERSE ] = 0;
        
        StringBuilder bases = new StringBuilder();
        StringBuilder quals = new StringBuilder();
        
        int i = -1;
        for( DataSpot spot : list )
        {
            ++i;
            if( null == spot )
                continue;
            i_spot.name = (String) key;            
            bases.append( spot.bases );
            quals.append( spot.quals );
            i_spot.read_length[ i ] = spot.bases.length();
            i_spot.read_name[ i ] = spot.bname;
        }

        i_spot.bases = bases.toString();
        i_spot.quals = quals.toString();
        i_spot.read_start[ IlluminaSpot.REVERSE ] = i_spot.read_length[ IlluminaSpot.FORWARD ];
        
        return i_spot;
    }

    
    @Override
    public boolean 
    isCollected( List<DataSpot> list )
    {
        return null != list.get( 0 ) 
               && null != list.get( 1 );
    }


    @Override
    public IlluminaSpot 
    handleErrors( final Object key, List<DataSpot> list ) throws DataConsumerException
    {   
        return assemble( key, list );
    }
}
