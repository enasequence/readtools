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
package uk.ac.ebi.ena.readtools.sampler;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.zip.GZIPInputStream;

import uk.ac.ebi.ena.readtools.loader.common.QualityNormalizer;
import uk.ac.ebi.ena.readtools.loader.common.feeder.AbstractDataProducer;
import uk.ac.ebi.ena.readtools.loader.common.feeder.DataProducerException;
import uk.ac.ebi.ena.readtools.loader.fastq.DataSpot;
import uk.ac.ebi.ena.readtools.loader.fastq.DataSpot.DataSpotParams;
import uk.ac.ebi.ena.readtools.sampler.intervals.QualInterval;
import uk.ac.ebi.ena.readtools.sampler.intervals.QualScoring;

public class 
FastqFileReader 
{
    SamplingDataEater<DataSpot, Object>eater = new SamplingDataEater<DataSpot, Object>();
    File file;

    
    public
    FastqFileReader( File file )
    {
        this.file = file;
    }
    
    
    private boolean 
    read( InputStream             is, 
          final QualityNormalizer normalizer,
          final String stream_marker ) throws SecurityException, DataProducerException, NoSuchMethodException, IOException
    {
        AbstractDataProducer<DataSpot> df = new AbstractDataProducer<DataSpot>( is, DataSpot.class )
        {
            final DataSpotParams params = DataSpot.defaultParams(); 
            
            @Override
            protected DataSpot
            newProducible()
            {
                return new DataSpot( normalizer, stream_marker, params );
            }
        };
        
        df.setConsumer( eater );
        df.run();
        return df.isOk();
    }
    
 /*   
    private boolean 
    read( String                   resource, 
          final  QualityNormalizer normalizer ) throws Exception
    {
        InputStream is = getClass().getResourceAsStream( resource );
        try
        {
            return read( is, normalizer );
        } finally
        {
            is.close();
        }
    }
*/    
    
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
        return read( file, QualityNormalizer.NONE );
    }
    
    
    public void
    printDataRow()
    {
        int lbound = eater.getLBoundary();
        System.out.println( String.format( "Reads: %d", eater.getFeedCount() ) );
        long[] interval = eater.getDataInterval();
        
        for( int i = 0; i < interval.length; ++ i )
            System.out.println( String.format( "%3d[%s]: %d", lbound + i, (char)( lbound + i ), interval[ i ] ) );
        
        
    }
    
    
    double
    getPhredP( int q )
    {
        return Math.pow( 10, - q / 10d );
    }


    double
    getLogOddsP( int q )
    {
        return Math.pow( 10, - q / 10d ) / ( 1 + Math.pow( 10, - q / 10d ) );
    }

    
    public void
    printDataRow3()
    {
        int lbound = eater.getLBoundary();
        long[] interval = eater.getBaseDataInterval();
        int peak_pos = -1;

        
        
        long peak_value = -1;
        for( int i = 0; i < interval.length; ++ i )
        {        
            if( interval[ i ] > peak_value )
            {
                peak_pos = i;
                peak_value = interval[ i ]; 
            }
        }


        System.out.println( String.format( "Bases: %d, peak: %3d[%s]", 
                                           eater.getBaseCount(),
                                           lbound + peak_pos,
                                           (char)( lbound + peak_pos ) ) );
        
        for( int i = 0; i < interval.length; ++ i )
        {
            if( interval[ i ] <= 0 )
                continue;
            
            System.out.print( String.format( "%3d[%s] %9d\t: ", 
                                             lbound + i, 
                                             (char)( lbound + i ),
                                             interval[ i ] ) );

            int to = (int) Math.log( Math.floor( (double)interval[ i ] /*/ (double)row_count */) );
            for( int j = 0; j < to; ++j )
                System.out.print( "*" );
            System.out.println();
        }
        
        
        long[][] tinterval = eater.getBaseTransitionDataInterval();
        List<Pair> sorted = new ArrayList<Pair>();
        
        for( int i = 0; i < tinterval.length; ++ i )
        {
            for( int k = 0; k < tinterval[ i ].length; ++ k )
            {
                if( tinterval[ i ][ k ] <= 0 )
                    continue;
                
                System.out.print( String.format( "%3d->%3d[%s->%s] %9d\t: ", 
                                                 lbound + i,
                                                 lbound + k, 
                                                 (char)( lbound + i ),
                                                 (char)( lbound + k ),
                                                 tinterval[ i ][ k ] ) );
    
                int to = (int) Math.log( Math.floor( (double)tinterval[ i ][ k ] /*/ (double)row_count */) );
                for( int j = 0; j < to; ++j )
                    System.out.print( "*" );
                System.out.println();
                sorted.add( new Pair( tinterval[ i ][ k ], String.format( "%3d->%3d[%s->%s]", 
                                                                          lbound + i,
                                                                          lbound + k, 
                                                                          (char)( lbound + i ),
                                                                          (char)( lbound + k ) ) ) );
            }
        }
        
        Collections.sort( sorted );
        for( Pair p : sorted )
            System.out.printf( "%s %9d\n",
                               p.value,
                               p.count );


    }



    public void
    printDataRow2( int offset )
    {
        int lbound = eater.getLBoundary();
        long row_count = eater.getFeedCount(); 
        
        long[] interval = eater.getDataInterval();
        
        int from_interval = 0;
        int to_interval   = interval.length;
        int peak_pos = -1;
        
        for( int i = 0; i < interval.length; ++ i )
        {        
            if( interval[ i ] > 0 )
            {
                from_interval = i;
                break;
            }
        }
        
        
        for( int i = interval.length -1; i >= 0; -- i )        
        {
            if( interval[ i ] > 0 )
            {
                to_interval = i;
                break;
            }
        }
        
        
        long peak_value = -1;
        for( int i = 0; i < interval.length; ++ i )
        {        
            if( interval[ i ] > peak_value )
            {
                peak_pos = i;
                peak_value = interval[ i ]; 
            }
        }


        System.out.println( String.format( "Reads: %d, lbound: %3d[%s], hbound: %3d[%s], peak: %3d[%s]", 
                                           eater.getFeedCount(),
                                           lbound + from_interval, 
                                           (char)( lbound + from_interval ),
                                           lbound + to_interval,
                                           (char)( lbound + to_interval ),
                                           lbound + peak_pos,
                                           (char)( lbound + peak_pos ) ) );
        
        List<Pair> sorted = new ArrayList<Pair>();
        
        for( int i = from_interval; i <= to_interval; ++ i )
        {
            System.out.print( String.format( "%3d[%s]%3d: %9d: \t%f \t%f\t: ", 
                                             lbound + i, 
                                             (char)( lbound + i ),
                                             lbound + i - offset, 
                                             interval[ i ],
                                             getPhredP( lbound + i - offset ),
                                             getLogOddsP( lbound + i - offset ) ) );
            int to = (int) Math.log( Math.floor( (double)interval[ i ] /*/ (double)row_count */) );
            for( int j = 0; j < to; ++j )
                System.out.print( "*" );
            System.out.println();
            
            sorted.add( new Pair( interval[ i ], String.format( "%3d", lbound + i ) ) );
        }
        
        Collections.sort( sorted );
        int index = 0;
        for( Pair p : sorted )
            System.out.printf( "case %s: return %3d; //%s[%s]: %9d\n",
                               p.value,
                               index++,
                               
                               p.value,
                               Character.valueOf( (char) Integer.parseInt( p.value.trim() ) ), 
                               p.count );

        
        
        
        long[][] tinterval = eater.getQualTransitionDataInterval();
        List<Pair> qsorted = new ArrayList<Pair>();
        long tr_count = 0;
        for( int i = 0; i < tinterval.length; ++ i )
        {
            for( int k = 0; k < tinterval[ i ].length; ++ k )
            {
                if( tinterval[ i ][ k ] <= 0 )
                    continue;
                tr_count += tinterval[ i ][ k ]; 
                System.out.print( String.format( "%3d->%3d[%s->%s] %9d\t: ", 
                                                 lbound + i,
                                                 lbound + k, 
                                                 (char)( lbound + i ),
                                                 (char)( lbound + k ),
                                                 tinterval[ i ][ k ] ) );
    
                int to = (int) Math.log( Math.floor( (double)tinterval[ i ][ k ] /*/ (double)row_count */) );
                for( int j = 0; j < to; ++j )
                    System.out.print( "*" );
                System.out.println();
                qsorted.add( new Pair( tinterval[ i ][ k ], String.format( "%3d->%3d[%s->%s]", 
                                                                          lbound + i,
                                                                          lbound + k, 
                                                                          (char)( lbound + i ),
                                                                          (char)( lbound + k ) ) ) );
            }
        }
        
        Collections.sort( qsorted );
        for( Pair p : qsorted )
        {
            System.out.printf( "%s %9d\t%12f\t",
                               p.value,
                               p.count,
                               (double)p.count / (double)tr_count );

            int to = (int) Math.log( Math.floor( (double)p.count /*/ (double)row_count */) );
            for( int j = 0; j < to; ++j )
                System.out.print( "*" );
            System.out.println();
        }


        
    }


    public QualInterval
    printIntervals()
    {
        Map<QualInterval, Double> list = QualInterval.forRow( eater.getDataInterval(), eater.getLBoundary() );
        double weight = 0;
        QualInterval result = null;
        for( Iterator<Entry<QualInterval, Double>> i = list.entrySet().iterator(); i.hasNext(); )
        {
            Entry<QualInterval, Double> e = i.next();
            System.out.println( e );
            if( weight < e.getValue() )
            {
                result = e.getKey();
                weight = e.getValue();
            }
        }
        
        return result;
    }
    
    
    public QualScoring
    printScoring( int offset )
    {
        Map<QualScoring, Double> map = QualScoring.forRow( eater.getDataInterval(), offset - eater.getLBoundary() );
        double weight = 126;
        QualScoring result = null;
        for( Iterator<Entry<QualScoring, Double>> i = map.entrySet().iterator(); i.hasNext(); )
        {
            Entry<QualScoring, Double> e = i.next();
            System.out.println( e );
            if( weight > e.getValue() )
            {
                result = e.getKey();
                weight = e.getValue();
            }
        }
        
        return result;
    }


    class
    Pair implements Comparable<Pair>
    {
        public String value;
        public long   count;
        
        
        Pair( long count, String value )
        {
            this.value = value;
            this.count = count;
        }

        public int 
        compareTo( Pair o )
        {
            return (int) -(this.count - o.count);
        }
    }
    

}
