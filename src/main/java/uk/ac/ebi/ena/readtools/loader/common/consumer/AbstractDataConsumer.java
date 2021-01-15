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
package uk.ac.ebi.ena.readtools.loader.common.consumer;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;


public abstract class
AbstractDataConsumer<T1 extends Spot, T2 extends Spot> implements DataConsumer<T1, T2>
{
    protected Map<String, List<T1>> spots = null;
    protected DataConsumer<T2, ?> dataConsumer;
    
    private long log_time =  System.currentTimeMillis();
    private long log_interval = 60 * 1000;
    private long assembled = 0;
    private long ate     = 0;
    protected boolean verbose = false;
    private volatile boolean is_ok = true;
    
    
    @Override
    public boolean 
    isOk()
    {
        return null == dataConsumer ? is_ok : is_ok && dataConsumer.isOk();
    }
    
    
    public AbstractDataConsumer<T1, T2>
    setVerbose( boolean verbose )
    {
        this.verbose = verbose;
        return this;
    }
    
    
    public AbstractDataConsumer()
    {
        this( 1024 * 1024 ); 
    }


    protected AbstractDataConsumer(int map_size )
    {
        spots = new HashMap<>( map_size );
    }
    
    
    @Override
    public void
    setConsumer(DataConsumer<T2, ? extends Spot> dataConsumer)
    {
        if( !dataConsumer.equals( this ) )
            this.dataConsumer = dataConsumer;
    }

    
    public abstract String
    getKey( T1 spot ) throws DataConsumerException;
    
    
    public abstract T2
    assemble(final String key, List<T1> list ) throws DataConsumerException;
    
    
    public void 
    append( List<T1> list, T1 spot ) throws DataConsumerException
    {
        list.add( spot );
    }
    
    
    public List<T1>
    newListBucket()
    {
        return Collections.synchronizedList( new ArrayList<T1>() );
    }
    
    
    public abstract boolean
    isCollected( List<T1> list );
    
    
    public synchronized void
    cascadeErrors() throws DataConsumerException
    {
        for( Entry<String, List<T1>> entry : spots.entrySet() )
        {
            if( null != dataConsumer)
                dataConsumer.consume( handleErrors( entry.getKey(), entry.getValue() ) );
            else 
                System.out.println( "<?> " + handleErrors( entry.getKey(), entry.getValue() ) );
        }
        
        if( null != dataConsumer)
            dataConsumer.cascadeErrors();
    }
    
    
    public abstract T2
    handleErrors(final String key, List<T1> list ) throws DataConsumerException;
    
    
    public void
    consume(T1 spot ) throws DataConsumerException
    {
        //System.out.println( spot );
        
        String key = getKey( spot );
        List<T1> list = null;
        
        synchronized(spots)
        {
            if( !spots.containsKey( key )  )
            {
                list = newListBucket();
                spots.put( key, list );
            } else
            {
                list = spots.get( key );
            }

            ate ++;
            append( list, spot );
        
            if( isCollected( list ) )
            {
                T2 assembly = assemble( key, list );
                assembled ++;
 //               synchronized( spots )
                {
                    spots.remove( key );
                }
                
                if( null != dataConsumer)
                    dataConsumer.consume( assembly );
                else
                    System.out.println( assembly );
            }
        }
        
        
        if( verbose )
        {
            long time = System.currentTimeMillis();
        
            if( time > log_time )
            {
                log_time = time + log_interval;
                System.out.println( String.format( "Ate: %d,\tAssembled: %d,\tDiff: %d", ate, assembled, ate - ( assembled << 1 ) ) );
            }
        }

    }
}
