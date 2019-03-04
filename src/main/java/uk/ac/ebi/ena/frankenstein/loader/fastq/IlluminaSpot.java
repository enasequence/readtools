package uk.ac.ebi.ena.frankenstein.loader.fastq;



public class 
IlluminaSpot
{
    //Must be private. DO NOT USE!
    private 
    IlluminaSpot()
    {
        ;
    }
    
    
    public static final int FORWARD = 0;
    public static final int REVERSE = 1;

    //NCBI coords start from 0!
    public int[] read_start;
    public int[] read_length;
    //TODO byte arrays?
    public String bases = "";
    public String quals = "";
    public String name  = "";
    
    
    public static IlluminaSpot 
    initPaired()
    {
        IlluminaSpot result = new IlluminaSpot();
        result.read_start = new int[] { -1, -1 };
        result.read_length = new int[] { -1, -1 };
        return result;
    }
    
    
    public static IlluminaSpot 
    initSingle()
    {
        IlluminaSpot result = new IlluminaSpot();
        result.read_start = new int[] { -1 };
        result.read_length = new int[] { -1 };
        return result;
    }
    
    
    public String
    toString()
    {
        StringBuilder result = new StringBuilder()
                               .append( name )
                               .append( '\n' )
                               .append( read_start[ FORWARD ] )
                               .append( ":" )
                               .append( read_length[ FORWARD ] );
        
        if( 2 == read_start.length )
            result.append( ", " )
                  .append( read_start[ REVERSE ] )
                  .append( ":" )
                  .append( read_length[ REVERSE ] );
                   
        result.append( '\n' )
              .append( bases )
              .append( '\n' )
              .append( quals );
        
        return result.toString();
    }
}
