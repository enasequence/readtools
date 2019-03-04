package uk.ac.ebi.ena.frankenstein.loader.cg.reads;

import java.util.Arrays;

public class 
CGSpot
{
    public enum 
    DNB
    {
        r1( 0,   5, DNBArm.LEFT ),
        r2( 5,  10, DNBArm.LEFT ),
        r3( 15, 10, DNBArm.LEFT ),
        r4( 25, 10, DNBArm.LEFT ),
        r5( 35, 10, DNBArm.RIGHT ),
        r6( 45, 10, DNBArm.RIGHT ),
        r7( 55, 10, DNBArm.RIGHT ),
        r8( 65,  5, DNBArm.RIGHT );

        DNB( int pos, int length, DNBArm half )
        {
            this.pos      = pos;
            this.length   = length;
            this.half     = half;
        }
        
        int pos    = -1;
        int length = -1;
        DNBArm half;
        
        public int
        getPos()
        {
            return pos;
        }
        
        
        public int
        getLength()
        {
            return length;
        }
        
        public DNBArm
        getHalf()
        {
           return half; 
        }
    }
    
    
    public static class
    CGMapping
    {
        String  chromosome;
        Integer offset;
        Integer weight;
        DNBArm  side;
        Strand  strand;
        Integer gaps[] = new Integer[ 3 ];
        Integer mateRec;
        
        public Integer[]
        getGaps()
        {
            return gaps;
        }
        
        public String
        getChromosome()
        {
            return chromosome;
        }
        
        public Integer
        getOffset()
        {
            return offset;
        }
        
        public Integer
        getWeight()
        {
            return weight;
        }
        
        public Strand
        getStrand()
        {
            return strand;
        }
        
        public Integer
        getMateRec()
        {
            return mateRec;
        }
        
        public DNBArm
        getSide()
        {
            return side;
        }
        
        public String
        toString()
        {
            return String.format( "%s:%d:%d:%s:%s:%s:%d", 
                                   chromosome, 
                                   offset, 
                                   weight,
                                   side,
                                   strand, 
                                   Arrays.asList( gaps ),
                                   mateRec );
        }
    }
    
    
    public enum
    DNBArm
    {
        LEFT,
        RIGHT
    }   
        

    public enum
    Strand
    {
        FORWARD,
        REVERSE
    }

    
    String    reads;
    String    scores;
    CGMapping mappings[];
    

    public String
    toString()
    {
        return String.format( "%s %s %s",
                              reads,
                              scores,
                              null == mappings ? "none" : Arrays.asList( mappings ) );
    }
    
}
