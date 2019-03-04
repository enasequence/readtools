package uk.ac.ebi.ena.sampler;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.zip.GZIPInputStream;

import uk.ac.ebi.ena.frankenstein.loader.common.QualityNormalizer;
import uk.ac.ebi.ena.frankenstein.loader.common.feeder.AbstractDataFeeder;
import uk.ac.ebi.ena.frankenstein.loader.common.feeder.DataFeederException;
import uk.ac.ebi.ena.frankenstein.loader.fastq.DataSpot;
import uk.ac.ebi.ena.frankenstein.loader.fastq.DataSpot.DataSpotParams;

public class 
SplitReader 
{
    StripesDataEater<DataSpot, Object>eater;
    File file;

    
    public
    SplitReader( File file, File folder ) throws IOException
    {
        this.file = file;
        eater = new StripesDataEater<DataSpot, Object>( folder );
    }
    
    
    private boolean 
    read( InputStream             is, 
          final QualityNormalizer normalizer,
          final String stream_marker ) throws SecurityException, DataFeederException, NoSuchMethodException, IOException
    {
        AbstractDataFeeder<DataSpot> df = new AbstractDataFeeder<DataSpot>( is, DataSpot.class ) 
        {
            final DataSpotParams params = DataSpot.defaultParams();
            
            @Override
            protected DataSpot newFeedable()
            {
                return new DataSpot( normalizer, stream_marker, params );
            }
        };
        
        df.setEater( eater );
        df.run();
        return df.isOk();
    }
    

    
    private boolean 
    read( File                    file, 
          final QualityNormalizer normalizer ) throws Exception
    {
        
        InputStream is = new BufferedInputStream( new FileInputStream( file ) );
        if( file.getName().endsWith( ".gz" ) )
            is = new GZIPInputStream( is );
        try
        {
            return read( is, normalizer, file.getName() );
        } finally
        {
            is.close();
        }
    }
    
    
    public boolean 
    process() throws Exception
    {
        boolean result = read( file, QualityNormalizer.NONE );
        eater.cascadeErrors();
        return result;
    }
    
    
    public long 
    getProcessedSpotsCount()
    {
        return eater.getFeedCount();
    }
    
  }
