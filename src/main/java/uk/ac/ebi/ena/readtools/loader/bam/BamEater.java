package uk.ac.ebi.ena.readtools.loader.bam;

import java.io.File;
import java.util.List;

import htsjdk.samtools.SAMRecord;
import uk.ac.ebi.ena.readtools.loader.common.eater.AbstractPagedDataEater;
import uk.ac.ebi.ena.readtools.loader.common.eater.DataEater;
import uk.ac.ebi.ena.readtools.loader.common.eater.DataEaterException;
import uk.ac.ebi.ena.readtools.loader.fastq.IlluminaSpot;

public class 
BamEater extends AbstractPagedDataEater<BamSpot, IlluminaSpot> 
         implements DataEater<BamSpot, IlluminaSpot>
{
    private static final Object BAM_STAR = "*"; 
    private boolean paired_fetched = false;
    private boolean allow_paired = false;
    private boolean allow_empty  = false;
    private int     errors = 0;
    private boolean verbose = true;
    

    public
    BamEater( File tmp_root, boolean allow_paired, boolean allow_empty, int spill_page_size )
    {
        super( tmp_root, spill_page_size );
        this.allow_paired = allow_paired;
        this.allow_empty = allow_empty;
    }
    
    
    @Override
    public void
    eat( BamSpot object ) throws DataEaterException
    {
        //do not load supplementary reads
        if( object.spot.getSupplementaryAlignmentFlag() )
            return;

        if( object.spot.getSecondaryAlignmentFlag() )
            return;
        
        if( object.spot.getDuplicateFlag() )
            return;
        
        //do not process empty reads
        if( allow_empty 
        	&& ( ( 0 == object.spot.getReadString().length() && 0 == object.spot.getBaseQualityString().length() ) 
        	     || ( 1 == object.spot.getReadString().length() && 1 == object.spot.getBaseQualityString().length() 
        	          && object.spot.getReadString().equals( BAM_STAR ) && object.spot.getReadString().equals( BAM_STAR ) ) ) )
        return;
        
        
        // allows only either paired or non-paired records
        if( !allow_paired && object.spot.getReadPairedFlag() )
            throw new UnsupportedOperationException(); //TODO: MetaDataException()3
        
        if( allow_paired )
            paired_fetched |= object.spot.getReadPairedFlag();
  
        super.eat( object );
    }
    
    
    @Override
    public Object 
    getKey( BamSpot object )
    {
        return object.spot.getReadName();
    }

    
    @Override
    public IlluminaSpot 
    assemble( final Object key, List<BamSpot> list )
    {
        IlluminaSpot s = null;
        if( !allow_paired )
        {
            s = IlluminaSpot.initSingle();
            BamSpot record = list.get( 0 );
            s.bases = record.spot.getReadString();
            s.quals = record.spot.getBaseQualityString();
            s.name  = record.spot.getReadName();
            s.read_length[ IlluminaSpot.FORWARD ] = s.bases.length();
            s.read_start[ IlluminaSpot.FORWARD ] = 0;
            
            if( record.spot.getReadNegativeStrandFlag() )
            {
                s.bases = new StringBuilder( new String( complement( s.bases.getBytes() ) ) ).reverse().toString();
                s.quals = new StringBuilder( s.quals ).reverse().toString();
            }
            
        } else 
        {
            // normalize
            if( list.size() == 1 )
                list.add( null );
            
            if( ( list.get( 0 ).spot.getReadPairedFlag() && !list.get( 0 ).spot.getFirstOfPairFlag() ) 
                || ( null != list.get( 1 ) && list.get( 1 ).spot.getFirstOfPairFlag() ) )
            {
                //swap
                BamSpot first = list.get( 0 );
                list.set( 0, list.get( 1 ) );
                list.set( 1, first );
            }
                
            s = IlluminaSpot.initPaired();
            s.read_start[ IlluminaSpot.FORWARD ] = 0;
            s.read_start[ IlluminaSpot.REVERSE ] = 0;

            s.read_length[ IlluminaSpot.FORWARD ] = 0;
            s.read_length[ IlluminaSpot.REVERSE ] = 0;

            
            StringBuilder bases = new StringBuilder();
            StringBuilder quals = new StringBuilder();
            
            int i = -1;
            for( BamSpot spot : list )
            {
                ++i;
                if( null == spot )
                    continue;
                
                //printSpot( spot, i + 1 );
                
                bases.append( spot.spot.getReadNegativeStrandFlag() ? new StringBuilder( new String( complement( spot.spot.getReadString().getBytes() ) ) ).reverse().toString() : spot.spot.getReadString() );
                quals.append( spot.spot.getReadNegativeStrandFlag() ? new StringBuilder( spot.spot.getBaseQualityString() ).reverse().toString() : spot.spot.getBaseQualityString() );
                
                s.name = (String)getKey( spot );
                s.read_length[ i ] = spot.spot.getReadString().length();
            }

            s.bases = bases.toString();
            s.quals = quals.toString();
            
            s.read_start[ IlluminaSpot.REVERSE ] = s.read_length[ IlluminaSpot.FORWARD ];
        } 
        
        if( null != s.name ) 
            s.name = s.name.replaceAll( " ", "_" );
        
        return s;
    }
    

    private static final byte a='a', c='c', g='g', t='t', A='A', C='C', G='G', T='T';


    public static byte 
    complement_table( final byte b ) 
    {
        switch ( b ) 
        {
            case a: return t;
            case c: return g;
            case g: return c;
            case t: return a;
            case A: return T;
            case C: return G;
            case G: return C;
            case T: return A;
            default: return b;
        }
    }

    public static byte[] 
    complement( final byte[] bases ) 
    {
        for( int i = 0; i < bases.length; ++i )
            bases[ i ] = complement_table( bases[ i ] );
        
        
        return bases;
    }

    
    private void
    printSpot( SAMRecord spot, int index )
    {
        int add = 36 - spot.getReadString().length();
        String b = "";
        String q = ""; 
//        if( add > 0 )
//        {
//            char[] buff = new char[ add ];
//            Arrays.fill( buff, 'A' ); 
//            b = new String( buff );
//            Arrays.fill( buff, '|' ); 
//            q = new String( buff );
//        } else if( add  < 0 )
//            throw new RuntimeException( "" + spot.getReadString().length() ); 
        String name = spot.getReadName().substring( 0, spot.getReadName().lastIndexOf( ':' ) );
        System.out.println( "@" + name + "#0/" + index );
        System.out.println( ( spot.getReadNegativeStrandFlag() ? new StringBuilder( spot.getReadString() ).reverse().toString() : spot.getReadString() ) + b );
        System.out.println( "+" + name + "#0/" + index );
        System.out.println( ( spot.getReadNegativeStrandFlag() ? new StringBuilder( spot.getBaseQualityString() ).reverse().toString() : spot.getBaseQualityString() ) + q );
        System.out.flush();
    }
    
    
    @Override
    public boolean 
    isCollected( List<BamSpot> list )
    {
        int size = list.size();
        if( allow_paired )
        {
                if( size == 2 )
                    return true;
        } else if( size == 1 )
            return true;
        return false;
    }

    
    @Override
    public IlluminaSpot 
    handleErrors( final Object key, List<BamSpot> list )
    {
        ++errors;
        return assemble( key, list );
    }

    
    public synchronized void
    cascadeErrors() throws DataEaterException
    {
        super.cascadeErrors();
        
        if( verbose )
            System.out.println( "Cascaded errors: " + errors );
    }
}
