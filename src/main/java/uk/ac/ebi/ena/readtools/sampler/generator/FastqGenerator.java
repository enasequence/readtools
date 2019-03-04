package uk.ac.ebi.ena.readtools.sampler.generator;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;

public class 
FastqGenerator
{
    public enum
    Scoring
    {
        PHRED(),
        LOGODDS();
        
        private char
        getPhredQ( double p )
        {
            return (char)( -10 * Math.log10( p ) );
        }
        
        
        private char 
        getLoggodsQ( double p )
        {
            return (char)( -10 * Math.log10( p / ( 1 - p ) ) );
        }

        
        char
        getQ( double p )
        {
            switch( this )
            {
            case PHRED:
                return getPhredQ( p );
                
            case LOGODDS:
                return getLoggodsQ( p );
                
            default:
                throw new RuntimeException( "Not Implemented" );
            }
        }
    }
    
    
    
    protected int lbound;
    protected int hbound;
    protected Scoring qual_scoring;
    protected int offset;
    protected int size;
    
    
    public static void
    main( String[] args )
    {
        Params params = new Params();
        JCommander jc = new JCommander( params );
        try
        {
            jc.parse( args );
        } catch ( Exception e )
        {
            jc.usage();
            System.exit( 1 );
        }
        
        
        FastqGenerator fq = new FastqGenerator();
        fq.configure( params );
        
        fq.produce( params.count );
    }
    
    
    void
    configure( Params params )
    {
        lbound = params.lboundary;
        hbound = params.hboundary;
        qual_scoring = Scoring.valueOf( params.quality_scoring.toUpperCase() );
        offset = params.quality_offset;
        size   = params.sequence_size;
    }
    
    
    protected void
    produce( int count )
    {
        
        for( int i = 0; i < count; ++ i )
        {
            String readname = getReadname( i ); 
            emit( readname, getRead( i ), getQuals( i ) );
        }
    }
    
    
    protected String
    getReadname( int number )
    {
        return String.format( "FASTQ_GENERATOR:%09d/1", number );
    }
    

    
    protected String
    getQuals( int number )
    {
        char value[] = new char[ size ];
        for( int i = 0; i < value.length; ++ i )
        { 
            int q = 6;
            do
            {
                q = qual_scoring.getQ( Math.random() );
            } while( q < lbound || q > hbound );
            
            int v = (int) ( q + offset );
            if( v < '!' || v > '~')
                throw new RuntimeException( "range exceeded: "  + v );
            value[ i ] = (char)v;
        }

        return new String( value );
    }

    
    static char[] BASES = { 'A', 'C', 'G', 'T' };
    
    protected String
    getRead( int number )
    {
        char value[] = new char[ size ];
        for( int i = 0; i < value.length; ++ i )
            value[ i ] = (char) BASES[ (int)( BASES.length * Math.random() ) ];
        
        return new String( value );
    }
    
    
    protected void
    emit( String readname, 
          String reads, 
          String quals )
    {
        System.out.printf( "@%s\n%s\n+%s\n%s\n", readname, reads, readname, quals );
    }
    
    
    
    static class
    Params
    {
        
        @Parameter( names = { "--quality-scoring", "-qs" }, description = "Type of quality scoring system", required = true )
        String quality_scoring;
        
        @Parameter( names = { "--quality-offset", "-qo" }, description = "Quality offset", required = true )
        int quality_offset;
        
        @Parameter( names = { "--low-boundary", "-lb" }, description = "Low boundary of quality interval in output", required = true )
        int lboundary;
        
        @Parameter( names = { "--high-boundary", "-hb" }, description = "High boundary of quality interval in output", required = true )
        int hboundary;
        
        @Parameter( names = { "--count", "-c" }, description = "Number of records to be produced", required = true )
        int count;
        
        @Parameter( names = { "--sequence-size", "-ss" }, description = "Size of sequences", required = true )
        int sequence_size;
    }
}
