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
package uk.ac.ebi.ena.readtools.sampler.intervals;

import java.util.HashMap;
import java.util.Map;

public enum 
QualInterval
{
    S( 33, 126, 33, 73, 33 ),
    X( 59, 126, 59, 104, 64 ),
    I( 64, 126, 64, 104, 64 ),
    J( 66, 126, 67, 104, 64 ),
    L( 33, 126, 33, 74, 33 );
    
    
    
    final protected int lbound;
    final protected int hbound;

    final protected int typical_lbound;
    final protected int typical_hbound;
    final protected int offset;
    
    private 
    QualInterval( int lbound, int hbound, int typical_lbound, int typical_hbound, int offset ) 
    {
        this.lbound = lbound;
        this.hbound = hbound;

        this.typical_lbound = typical_lbound;
        this.typical_hbound = typical_hbound;
        this.offset = offset;
    }
    
    
    double 
    getWeight( long[] interval, int offset )
    {
        long tcount = 0;
        long icount = 0;
        for( int i = 0; i < interval.length; ++ i )
        {        
            long count = interval[ i ];
            tcount += count;
            int pos = i + offset;
            if( count > 0 
                && pos >= typical_lbound
                && pos <= typical_hbound )
            {
                icount += count;
            }
        }

        return tcount == 0 ? 0 : (double) icount / (double)tcount;
    }
    
    
    public static Map<QualInterval, Double>
    forRow( long[] interval, int offset )
    {
        int from_interval = 0;
        int to_interval   = interval.length;

        
        for( int i = 0; i < interval.length; ++ i )
        {        
            if( interval[ i ] > 0 )
            {
                from_interval = i;
                break;
            }
        }
        
        
        for( int i = interval.length -1; i >= 0; -- i )        
        {
            if( interval[ i ] > 0 )
            {
                to_interval = i;
                break;
            }
        }
        
        from_interval += offset;
        to_interval   += offset;
        
        Map<QualInterval, Double> map = new HashMap<QualInterval, Double>();
        for( QualInterval i : QualInterval.values() )
        {
            if( from_interval >= i.lbound 
                && to_interval <= i.hbound )
                map.put( i, i.getWeight( interval, offset ) );
        }   

        return map;
    }
    
 
    public String
    toString()
    {
        return String.format( "%s( %d, %d, offset: %d )", super.toString(), lbound, hbound, offset );
        
    }
    
    
    public int
    getOffset()
    {
        return offset;
    }
}
