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
package uk.ac.ebi.ena.readtools.loader.common;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;

import uk.ac.ebi.ena.readtools.loader.common.QualityNormalizer.QualityNormaizationException;

public enum 
QualityConverter
{
    PHRED,
    LOGODDS;
    
    
    public byte[]
    convertNormalized( QualityConverter to, byte[] what )
    {
        byte[] result = new byte[ what.length ];
        if( this.equals( PHRED ) && to.equals( LOGODDS ) )
        {
            for( int i = 0; i < what.length; ++ i )
                result[ i ] = phred2solexa( what[ i ] );
        }else if( this.equals( LOGODDS ) && to.equals( PHRED ) )
        {
            for( int i = 0; i < what.length; ++ i )
                result[ i ] = solexa2phred( what[ i ] );
        
        }else if( this.equals( PHRED ) && to.equals( PHRED ) 
                  || this.equals( LOGODDS ) && to.equals( LOGODDS ) )
        {
            System.arraycopy( what, 0, result, 0, what.length );
        }
    
        return result;
    }
    
    final static int SOLEXA_MIN = -5;
    final static int PHRED_MIN  = 0;
    
    byte
    phred2solexa( byte phred )
    {
        if( phred == PHRED_MIN )
            return SOLEXA_MIN;
        else
        {
            byte result = (byte)Math.round( 10 * Math.log10( Math.pow( 10, (double)phred / 10 ) - 1 ) );
            return SOLEXA_MIN > result ? SOLEXA_MIN : result;
        }
    }
    
    
    byte
    solexa2phred( byte logodds )
    {
        return (byte)Math.round( 10 * Math.log10( Math.pow( 10, (double)logodds / 10 ) + 1 ) );
    }




    public static void
    main( String[] args ) throws QualityNormaizationException
    {
        class Params
        {
            @Parameter( names={ "-name_from" } )
            String name_from;
            
            @Parameter( names={ "-name_to" } )
            String name_to;

            @Parameter( names={ "-from" } )
            String from;
            
            @Parameter( names={ "-to" } )
            String to;
            
            @Parameter( names={ "-line" } )
            String line;
            
            @Parameter( names={ "-h", "--help" } )
            boolean help;
        }
        
        
        Params p = new Params();
        JCommander jc = new JCommander( p );
        jc.parse( args );
        
        if( p.help )
        {
            jc.usage();
            System.exit( 0 );
        }
        
        System.out.println( "FROM: " + p.name_from + ", " + p.from  );
        System.out.println( "TO:   " + p.name_to   + ", " + p.to  );
        System.out.println( p.line );
/*
        byte[] phred = { 20, 19, 18, 17, 16, 
                15, 14, 13, 12, 11,
                10,  9,  8,  7,  6,
                 5,  4,  3,  2,  1,
                 0 };
*/
        
        byte[] normalized = QualityNormalizer.valueOf( p.name_from ).normalize( p.line );
        byte[] converted  = QualityConverter.valueOf( p.from ).convertNormalized( QualityConverter.valueOf( p.to ), normalized );
        byte[] bconverted = QualityConverter.valueOf( p.to ).convertNormalized( QualityConverter.valueOf( p.from ), converted );

        for( int i = 0; i < normalized.length; ++i )
            System.out.println( String.format( "%d\t%d\t%d", normalized[ i ], converted[ i ], bconverted[ i ] ) );
        
        System.out.println( QualityNormalizer.valueOf( p.name_to ).denormalize( QualityConverter.valueOf( p.from ).convertNormalized( QualityConverter.valueOf( p.to ),normalized ) ) );
        
    }

}
