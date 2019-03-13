package uk.ac.ebi.ena.readtools.loader.common.eater;



public class 
NullDataEater<T1> implements DataEater<T1, Void>
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
        //do nothing, this is null data eater
    }

    
    @Override public void 
    setEater( DataEater<Void, ?> dataEater )
    {
        throw new RuntimeException( "Not implemented" );
    }


    @Override public boolean 
    isOk()
    {
        return true;
    }

}
