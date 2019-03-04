package uk.ac.ebi.ena.frankenstein.loader.cg.reads;

import java.util.List;

import uk.ac.ebi.ena.frankenstein.loader.cg.reads.CGSpot.CGMapping;
import uk.ac.ebi.ena.frankenstein.loader.common.eater.AbstractDataEater;
import uk.ac.ebi.ena.frankenstein.loader.common.eater.DataEaterException;

public class 
CompleteGenomicsPairedEater extends AbstractDataEater<CompleteGenomicsBase, CGSpot>
{

    @Override
    public Object 
    getKey( CompleteGenomicsBase object ) throws DataEaterException
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
    assemble( final Object key, List<CompleteGenomicsBase> list ) throws DataEaterException
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
                throw new DataEaterException();
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
    handleErrors( final Object key, List<CompleteGenomicsBase> list ) throws DataEaterException
    {
        // TODO Auto-generated method stub
        return null;
    }

}
