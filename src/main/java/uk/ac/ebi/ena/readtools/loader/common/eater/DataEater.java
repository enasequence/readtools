package uk.ac.ebi.ena.readtools.loader.common.eater;

public interface 
DataEater<T1, T2>
{
    void    cascadeErrors() throws DataEaterException;
    void    eat( T1 object ) throws DataEaterException;
    void    setEater( DataEater<T2, ?> dataEater );
    boolean isOk();
}
