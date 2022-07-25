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
package uk.ac.ebi.ena.readtools.webin.cli.rawreads;

import com.google.common.hash.BloomFilter;
import com.google.common.hash.Funnels;

import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicLong;
import java.util.stream.Collectors;

public class 
BloomWrapper 
{
    private static final double falsePositiveProbability = 0.01;

    //Bloom bloom;
    private BloomFilter<String> bloom;
    private AtomicLong addCount = new AtomicLong();
    private AtomicLong possibleDuplicateCount = new AtomicLong();
    private int possibleDuplicatesRetainLimit;
    private Set<String> possibleDuplicates;
    
    
    public 
    BloomWrapper( long expectedReads )
    {
        this( expectedReads, 25 );
    }


    /**
     *
     * @param expectedReads Expected number of reads that will be added into this instance.
     * @param possibleDuplicatesRetainLimit Maximum number of possible duplicates to retain in memory for reporting and
     *                                      duplicate verification. Once the limit is reached, more possible
     *                                      duplicate reads names will be dropped.
     */
    BloomWrapper( long expectedReads, int possibleDuplicatesRetainLimit)
    {
        this.possibleDuplicatesRetainLimit = possibleDuplicatesRetainLimit;
        this.bloom = BloomFilter.create( Funnels.unencodedCharsFunnel(), expectedReads, falsePositiveProbability);
        possibleDuplicates = new HashSet<>( this.possibleDuplicatesRetainLimit);
    }


    public void
    add( String readName )
    {
        addCount.incrementAndGet();
        
        if( bloom.mightContain( readName ) )
        {
            possibleDuplicateCount.incrementAndGet();
            if( possibleDuplicates.size() < possibleDuplicatesRetainLimit)
            {
                possibleDuplicates.add( readName );
            }
        } else
        {
            bloom.put( readName );
        }
    }


    /**
     * @return Number of times read names were added into this instance.
     */
    public long
    getAddCount()
    {
        return addCount.get();
    }
    
    
    public boolean
    hasPossibleDuplicates()
    {
        return getPossibleDuplicateCount() > 0;
    }
    
    
    public Long
    getPossibleDuplicateCount()
    {
        return possibleDuplicateCount.get();
    }

    
    public Set<String>
    getPossibleDuplicates()
    {
        return possibleDuplicates;
    }
    
    
    public Map<String, Set<String>>
    findAllduplications( String[] read_names, int limit )
    {
        Map<String, Integer>      counts = new HashMap<>( limit );
        Map<String, Set<String>>  result = new HashMap<>( limit );
        
        if( !hasPossibleDuplicates() )
            return result;
        
        AtomicInteger index = new AtomicInteger();
        for( String read_name : read_names )
        {
            index.incrementAndGet();
            if( contains( read_name ) )
            {
                counts.put( read_name, counts.getOrDefault( read_name, 0 ) + 1 );
                Set<String> dlist = result.getOrDefault( read_name, new LinkedHashSet<>() );
                dlist.add( String.valueOf(index.get()));
                result.put( read_name, dlist );
            }
            //TODO
            if( result.size() >= limit )
                break;
        }
        //return result;
        return result.entrySet().stream().filter( e-> counts.get(e.getKey()) > 1 ).collect( Collectors.toMap(e -> e.getKey(), e -> e.getValue() ) );
    }
    
   
    public boolean 
    contains( String readName )
    {
        return bloom.mightContain( readName );
    }

    public BloomWrapper getCopy() {
        BloomWrapper res = new BloomWrapper(1);
        res.bloom = this.bloom.copy();
        res.addCount = new AtomicLong(this.addCount.get());
        res.possibleDuplicateCount = new AtomicLong(this.possibleDuplicateCount.get());
        res.possibleDuplicatesRetainLimit = this.possibleDuplicatesRetainLimit;
        res.possibleDuplicates = new HashSet<>(res.possibleDuplicatesRetainLimit);
        res.possibleDuplicates.addAll(this.possibleDuplicates);

        return res;
    }
}
