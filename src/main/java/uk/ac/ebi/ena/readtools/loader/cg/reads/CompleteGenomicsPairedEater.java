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
package uk.ac.ebi.ena.readtools.loader.cg.reads;

import java.util.List;

import uk.ac.ebi.ena.readtools.loader.cg.reads.CGSpot.CGMapping;
import uk.ac.ebi.ena.readtools.loader.common.eater.AbstractDataConsumer;
import uk.ac.ebi.ena.readtools.loader.common.eater.DataConsumerException;

public class 
CompleteGenomicsPairedEater extends AbstractDataConsumer<CompleteGenomicsBase, CGSpot>
{

    @Override
    public Object 
    getKey( CompleteGenomicsBase object ) throws DataConsumerException
    {
        return object.getKey();
    }

    
    protected void
    fillReads( CompleteGenomicsRead read, CGSpot spot )
    {
        spot.reads  = read.reads;
        spot.scores = read.scores;
    }
    

    protected void
    fillMaps( CompleteGenomicsMap3 map, CGSpot.CGMapping mappings[] )
    {
        CGSpot.CGMapping m = new CGMapping();
        m.chromosome = map.chromosome;
        m.weight  = map.weight;
        m.offset  = map.offset;
        m.side    = ( map.flags & CompleteGenomicsMap3.Flags.side ) == 0 ? CGSpot.DNBArm.LEFT 
                                                                         : CGSpot.DNBArm.RIGHT;
        m.strand  = ( map.flags & CompleteGenomicsMap3.Flags.strand ) == 0 ? CGSpot.Strand.FORWARD 
                                                                           : CGSpot.Strand.REVERSE;
        m.gaps    = map.gaps;
        m.mateRec = map.mateRec;
        mappings[ map.batch_index ] = m;
    }

    
    
    @Override
    public CGSpot 
    assemble( final Object key, List<CompleteGenomicsBase> list ) throws DataConsumerException
    {
        CGSpot spot = new CGSpot();
        spot.mappings = new CGSpot.CGMapping[ list.size() - 1 ];  
        for( CompleteGenomicsBase l : list )
        {
            if( l instanceof CompleteGenomicsRead )
                fillReads( (CompleteGenomicsRead) l, spot );
            else if( l instanceof CompleteGenomicsMap3 )
                fillMaps( (CompleteGenomicsMap3) l, spot.mappings );
            else
                throw new DataConsumerException();
        }
        
        return spot;
    }

    @Override
    public boolean 
    isCollected( List<CompleteGenomicsBase> list )
    {
        boolean read = false;
        boolean maps = false;
        for( CompleteGenomicsBase l : list )
        {
            if( l instanceof CompleteGenomicsRead )
            {
                if( null == l.getKey() )
                    return true;
                
                read = true;
            } else if( l instanceof CompleteGenomicsMap3 )
            {
                CompleteGenomicsMap3 map = (CompleteGenomicsMap3) l;
                if( ( map.flags & CompleteGenomicsMap3.Flags.LastDNBRecord ) > 0 )
                    maps = true;
            } 
            
            if( read && maps )
                return true;
        }       
        
        return false;
    }

    @Override
    public CGSpot 
    handleErrors( final Object key, List<CompleteGenomicsBase> list ) throws DataConsumerException
    {
        // TODO Auto-generated method stub
        return null;
    }

}
