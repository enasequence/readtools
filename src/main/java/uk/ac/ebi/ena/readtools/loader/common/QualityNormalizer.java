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
package uk.ac.ebi.ena.readtools.loader.common;

/*
 * Normalises qualities read from file to NCBI VDB2 format
 */
public enum 
QualityNormalizer
{   /*
    S - Sanger        Phred+33,  raw reads typically (0, 40)
    X - Solexa        Solexa+64, raw reads typically (-5, 40)
    I - Illumina 1.3+ Phred+64,  raw reads typically (0, 40)
    J - Illumina 1.5+ Phred+64,  raw reads typically (3, 40)
    */
    
    SANGER( 33,  0, 93 /* 63 */ ),
    SOLEXA( 64, -5, 93 ),
    ILLUMINA_1_3( 64, 0, 93 /* 40 */ ),
    ILLUMINA_1_5( 64, 3, 93 ),

    //added to comply with sra-schema 
    X( 33, 0, 93 ),
    X_2( 64, 0, 93 ),
    
    // special type, leaves data as is for analytical usage.
    NONE( 33, 0, 93 );
    
    
    
    QualityNormalizer( int offset, int min, int max )
    {
        this.offset = (byte) offset;
        this.min    = (byte) min;
        this.max    = (byte) max;
    }
    

    byte offset;
    byte min;
    byte max;
    
    
    public byte[]
    normalize( String value ) throws QualityNormaizationException
    {
        if( min > max )
            throw new QualityNormaizationException( "min > max" );

        byte[] result = new byte[ value.length() ]; 
        for( int i = 0; i < result.length; ++i )
        {
            /* textual form with offset */
            byte q = (byte)( value.charAt( i ) - offset );
            if( q < min || q > max )
                throw new QualityNormaizationException( String.format(  "Quality value out of range, expected %s [%d, %d], offset %d, got %d (%s)", 
                                                                        this.name(), 
                                                                        this.min,
                                                                        this.max,
                                                                        this.offset,
                                                                        q,
                                                                        value.charAt( i ) ) );
            result[ i ] = q;
        }
        return result;
    }
    
    //
    // 
    public String
    denormalize( byte[] value ) throws QualityNormaizationException
    {
        if( min > max )
            throw new QualityNormaizationException( "min > max" );

        StringBuilder result = new StringBuilder();
        for( int i = 0; i < value.length; ++i )
        {
            /* textual to with offset */
            int q = value[ i ];
            if( q < min || q > max )
                throw new QualityNormaizationException( "Quality value out of range" );

            result.append( (char) ( q + offset ) );
        }
        return result.toString();
    }

    
    public class
    QualityNormaizationException extends Exception
    {
        private static final long serialVersionUID = 1L;

        public 
        QualityNormaizationException( String value )
        {
            super( value );
        }
    }
}
