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
package uk.ac.ebi.ena.readtools.loader.bam;

import java.io.IOException;
import java.io.InputStream;
import java.io.Serializable;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;

import uk.ac.ebi.ena.readtools.loader.common.feeder.FeedableData;

public class 
BamSpot implements Serializable
{
    private static final SAMFileHeader header = new SAMFileHeader();
    private static final long serialVersionUID = 5495870024344305947L;
    @FeedableData( method = "readSAMSpot" )
    MySAMRecord spot;
    
    
    static class 
    MySAMRecord
    {
        final private int    flags;
        final private String bases;
        final private String quals;
        final private String name; 
        

        private static final int READ_PAIRED_FLAG = 0x1;
        private static final int PROPER_PAIR_FLAG = 0x2;
        private static final int READ_UNMAPPED_FLAG = 0x4;
        private static final int MATE_UNMAPPED_FLAG = 0x8;
        private static final int READ_STRAND_FLAG = 0x10;
        private static final int MATE_STRAND_FLAG = 0x20;
        private static final int FIRST_OF_PAIR_FLAG = 0x40;
        private static final int SECOND_OF_PAIR_FLAG = 0x80;
        private static final int SECONDARY_ALIGNMENT_FLAG = 0x100;
        private static final int READ_FAILS_VENDOR_QUALITY_CHECK_FLAG = 0x200;
        private static final int DUPLICATE_READ_FLAG = 0x400;
        private static final int SUPPLEMENTARY_ALIGNMENT_FLAG = 0x800;
        
        
        public
        MySAMRecord( SAMRecord record )
        {
            this( record.getReadName(), record.getFlags(), record.getReadString(), record.getBaseQualityString() );  
        }
        
        
        public
        MySAMRecord( String name, int flags, String bases, String quals )
        {
            this.flags = flags;
            this.bases = bases;
            this.quals = quals;
            this.name  = name; 
        }
        
        public String
        getReadName()
        {
            return name;
        }
        
        
        public String
        getReadString()
        {
            return bases;
        }
        
        
        public String
        getBaseQualityString()
        {
            return quals;
        }
        
        
        public boolean
        getReadNegativeStrandFlag()
        {
            return 0 != ( flags & READ_STRAND_FLAG );
        }
        
        
        public boolean 
        getReadPairedFlag()
        {
            return 0 != ( flags & READ_PAIRED_FLAG );
        }
        
        
        public boolean
        getFirstOfPairFlag()
        {
            if( getReadPairedFlag() )
                return 0 != ( flags & FIRST_OF_PAIR_FLAG );
            throw new IllegalStateException( "Illegal method call for unpaired read!" );
        }
        
        
        public boolean
        getSupplementaryAlignmentFlag()
        {
            return 0 != ( flags & SUPPLEMENTARY_ALIGNMENT_FLAG );
        }
        
        
        public boolean
        getSecondaryAlignmentFlag()
        {
            return 0 != ( flags & SECONDARY_ALIGNMENT_FLAG );
        }
        
        
        public boolean
        getDuplicateFlag()
        {
            return 0 != ( flags & DUPLICATE_READ_FLAG );
        }
    }


    public
    BamSpot()
    {
        ;
    }
    
    
    public
    BamSpot( SAMRecord spot )
    {
        this.spot = new MySAMRecord( spot );
    }
    
    
    
    public void
    readSAMSpot( InputStream is )
    {
       if( true )
           throw new UnsupportedOperationException();
    }
    
    
    private void 
    writeObject( java.io.ObjectOutputStream out ) throws IOException
    {
        out.writeObject( spot.name );
        out.writeObject( spot.flags );
        out.writeObject( spot.bases );
        out.writeObject( spot.quals );
    }
    
    
    private void 
    readObject( java.io.ObjectInputStream in ) throws IOException, ClassNotFoundException
    {
        spot = new MySAMRecord( (String) in.readObject(), /* name  */ 
                                (Integer)in.readObject(), /* flags */
                                (String) in.readObject(), /* bases */
                                (String) in.readObject()  /* quals */
                               );
    }

  
}
