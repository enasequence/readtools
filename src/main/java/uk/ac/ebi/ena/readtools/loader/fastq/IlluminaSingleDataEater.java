package uk.ac.ebi.ena.readtools.loader.fastq;

import uk.ac.ebi.ena.readtools.loader.common.eater.DataEater;
import uk.ac.ebi.ena.readtools.loader.common.eater.DataEaterException;

public class 
IlluminaSingleDataEater implements DataEater<DataSpot, IlluminaSpot>
{
    DataEater<IlluminaSpot, ?> dataEater;
    boolean is_ok = true;
    
    @Override
    public void 
    cascadeErrors() throws DataEaterException
    {
        //empty, no errors here
        
        if( null != dataEater )
            dataEater.cascadeErrors();
    }

    
    @Override
    public void 
    eat( DataSpot spot ) throws DataEaterException
    {
        
        IlluminaSpot ispot = IlluminaSpot.initSingle();
        int slash_idx = spot.bname.lastIndexOf( '/' );
        ispot.name = slash_idx == -1 ? spot.bname 
                                     : spot.bname.substring( 0, slash_idx );
        ispot.bases = spot.bases;
        ispot.quals = spot.quals;
        ispot.read_start[ IlluminaSpot.FORWARD ] = 0;
        ispot.read_length[ IlluminaSpot.FORWARD ] = spot.bases.length();
        ispot.read_name[ IlluminaSpot.FORWARD ] = spot.bname;
        if( null != dataEater )
            dataEater.eat( ispot );
        else
            System.out.println( ispot );
    }


    @Override
    public void 
    setEater( DataEater<IlluminaSpot, ?> dataEater )
    {
        this.dataEater = dataEater;
    }


    @Override
    public boolean 
    isOk()
    {
        return null == dataEater ? is_ok : is_ok && dataEater.isOk();
    }
}
