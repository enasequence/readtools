package uk.ac.ebi.ena.readtools.loader.common.eater;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;


public abstract class 
AbstractDataEater<T1, T2> implements DataEater<T1, T2>
{
    protected Map<Object, List<T1>> objects = null; 
    protected DataEater<T2, ?> dataEater;
    
    private long log_time =  System.currentTimeMillis();
    private long log_interval = 60 * 1000;
    private long assembled = 0;
    private long ate     = 0;
    protected boolean verbose = false;
    private volatile boolean is_ok = true;
    
    
    protected void 
    setOk( boolean is_ok )
    {
        this.is_ok = is_ok;
    }
    
    
    @Override
    public boolean 
    isOk()
    {
        return null == dataEater ? is_ok : is_ok && dataEater.isOk();
    }
    
    
    public AbstractDataEater<T1, T2>
    setVerbose( boolean verbose )
    {
        this.verbose = verbose;
        return this;
    }
    
    
    public
    AbstractDataEater()
    {
        this( 1024 * 1024 ); 
    }


    protected
    AbstractDataEater( int map_size )    
    {
        objects = new HashMap<Object, List<T1>>( map_size );
    }
    
    
    @Override
    public void 
    setEater( DataEater<T2, ?> dataEater )
    {
        if( !dataEater.equals( this ) )
            this.dataEater = dataEater;
    }

    
    public abstract Object
    getKey( T1 object ) throws DataEaterException;
    
    
    public abstract T2
    assemble( final Object key, List<T1> list ) throws DataEaterException;
    
    
    public void 
    append( List<T1> list, T1 obj ) throws DataEaterException
    {
        list.add( obj );
    }
    
    
    public List<T1>
    newListBucket()
    {
        return Collections.synchronizedList( new ArrayList<T1>() );
    }
    
    
    public abstract boolean
    isCollected( List<T1> list );
    
    
    public synchronized void
    cascadeErrors() throws DataEaterException
    {
        for( Entry<Object, List<T1>> entry : objects.entrySet() )
        {
            if( null != dataEater )
                dataEater.eat( handleErrors( entry.getKey(), entry.getValue() ) );
            else 
                System.out.println( "<?> " + handleErrors( entry.getKey(), entry.getValue() ) );
        }
        
        if( null != dataEater )
            dataEater.cascadeErrors();
    }
    
    
    public abstract T2
    handleErrors( final Object key, List<T1> list ) throws DataEaterException;
    
    
    public void
    eat( T1 object ) throws DataEaterException
    {
        //System.out.println( object );
        
        Object key = getKey( object );
        List<T1> list = null;
        
        synchronized( objects )
        {
            if( !objects.containsKey( key )  )
            {
                list = newListBucket();
                objects.put( key, list );
            } else
            {
                list = objects.get( key );
            }

            ate ++;
            append( list, object );
        
            if( isCollected( list ) )
            {
                T2 assembly = assemble( key, list );
                assembled ++;
 //               synchronized( objects )
                {
                    objects.remove( key );
                }
                
                if( null != dataEater )
                    dataEater.eat( assembly );
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
