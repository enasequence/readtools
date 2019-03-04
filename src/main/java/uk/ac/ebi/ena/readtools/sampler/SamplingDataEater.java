package uk.ac.ebi.ena.readtools.sampler;


import uk.ac.ebi.ena.readtools.loader.common.QualityNormalizer;
import uk.ac.ebi.ena.readtools.loader.common.QualityNormalizer.QualityNormaizationException;
import uk.ac.ebi.ena.readtools.loader.common.eater.DataEater;
import uk.ac.ebi.ena.readtools.loader.common.eater.DataEaterException;
import uk.ac.ebi.ena.readtools.loader.fastq.DataSpot;

public class 
SamplingDataEater<T1 extends DataSpot, T2> implements DataEater<T1, T2>
{
    final protected int lbound = 33;  //!
    final protected int hbound = 126; //~;
    
    final protected long row[] = new long[ hbound - lbound + 1 ];
    final protected long rrow[] = new long[ hbound - lbound + 1 ];
    
    final protected long transition[][] = new long[ hbound - lbound + 1 ][ hbound - lbound + 1 ];
    final protected long rtransition[][] = new long[ hbound - lbound + 1 ][ hbound - lbound + 1 ];
    
    long  read_count;
    long  base_count;
    
    byte  q_prev = hbound - lbound - 1;
    byte  b_prev = hbound - 1;
    
    boolean lock;
    
    @Override
    public void 
    cascadeErrors() throws DataEaterException
    {
        throw new UnsupportedOperationException();
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

    
    @Override
    public void 
    eat( T1 object ) throws DataEaterException
    {
        try
        {
            for( byte q: QualityNormalizer.NONE.normalize( object.quals ) )
            {
                row[ q ] += 1;
                transition[ q_prev ][ q ] += 1;
                q_prev = q;
            }
            
            byte[] arr = object.bases.getBytes();
            for( int i = 0; i < arr.length; ++i )
            {
                byte b = arr[ i ];
                rrow[ b - lbound ] += 1;
                rtransition[ b_prev - lbound ][ b - lbound ] += 1;
                b_prev = b;
            }
            base_count += object.bases.length();
        } catch( QualityNormaizationException e )
        {
            System.out.printf( "read %d : %s, %s\n", read_count + 2, object, e.toString() );
        }
        
        ++read_count;
    }


    @Override
    public void 
    setEater( DataEater<T2, ?> dataEater )
    {
        throw new UnsupportedOperationException();
    }

    
    
    public long 
    getFeedCount()
    {
        return read_count;
    }
    

    public long 
    getBaseCount()
    {
        return base_count;
    }

    
    public long[]
    getBaseCounts()
    {
        return row;
    }
    
    
    public long[]
    getDataInterval()
    {
        return row;
    }
    
    
    public long[]
    getBaseDataInterval()
    {
        return rrow;
    }
                
    
    public long[][]
    getBaseTransitionDataInterval()
    {
        return rtransition;
    }
    
    public long[][]
    getQualTransitionDataInterval()
    {
        return transition;
    }

    
    public byte
    getLBoundary()
    {
        return lbound;
    }
    
    
    public byte
    getHBoundary()
    {
        return hbound;
    }


    @Override
    public boolean 
    isOk()
    {
        return true;
    }

}
