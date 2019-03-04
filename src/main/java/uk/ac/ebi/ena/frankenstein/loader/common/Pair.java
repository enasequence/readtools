package uk.ac.ebi.ena.frankenstein.loader.common;

import java.io.Serializable;

public class 
Pair<K, V> implements Serializable
{
    private static final long serialVersionUID = -4356524392295831193L;
    public K key;
    public V value;
    
    
    public 
    Pair()
    {
        this( null, null );
    }

    
    public 
    Pair( K key, V value )
    {
        this.key   = key;
        this.value = value;
    }
}
