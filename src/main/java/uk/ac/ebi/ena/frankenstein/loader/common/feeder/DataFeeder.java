package uk.ac.ebi.ena.frankenstein.loader.common.feeder;

import uk.ac.ebi.ena.frankenstein.loader.common.eater.DataEater;


public interface 
DataFeeder<T>
{
    public T feed() throws DataFeederEOFException, DataFeederException, DataFeederPanicException;
    public DataFeeder<T> setEater( DataEater<T, ?> eater );
    public boolean isOk();
    public Throwable getStoredException();
}
