package uk.ac.ebi.ena.frankenstein.loader.common.eater;



public class 
PrintDataEater<T1, T2> implements DataEater<T1, T2>
{
    volatile boolean lock;
    
    @Override
    public void 
    cascadeErrors() throws DataEaterException
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
    eat( T1 object ) throws DataEaterException
    {
        System.out.println( object ); 
    }

    
    @Override public void 
    setEater( DataEater<T2, ?> dataEater )
    {
        throw new UnsupportedOperationException();
    }


    @Override public boolean 
    isOk()
    {
        return true;
    }

}
