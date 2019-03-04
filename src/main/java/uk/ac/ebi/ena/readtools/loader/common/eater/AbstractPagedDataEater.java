package uk.ac.ebi.ena.readtools.loader.common.eater;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.EOFException;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.zip.Deflater;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import uk.ac.ebi.ena.readtools.loader.common.Pair;

public abstract class 
AbstractPagedDataEater<T1, T2> extends AbstractDataEater<T1, T2>
{
    private static final int OUTPUT_BUFFER_SIZE = 8192;
    private boolean use_spill = true;
    private List<File> files = new ArrayList<File>();    
    // note - must be at least n - 1; 
    final int spill_page_size;  //( Runtime.getRuntime().maxMemory() - Runtime.getRuntime().freeMemory() ) / 5120;
    final private File tmp_root;
    
    
    public
    AbstractPagedDataEater()
    {
        this( new File( "." ), 
              (int) ( ( Runtime.getRuntime().maxMemory() - Runtime.getRuntime().freeMemory() ) / 5120 ) );
    }

    
    public
    AbstractPagedDataEater( File tmp_root, 
                            int  spill_page_size )
    {
        super( spill_page_size );
        this.tmp_root = tmp_root;
        this.spill_page_size = spill_page_size;
    }
    
    
    private File
    getTempFile() throws IOException
    {
        String prefix = String.format( "TREAD_%d_TIME_%d_FILE_", Thread.currentThread().getId(), System.currentTimeMillis() );
        String suffix = String.format( "_PAGE_%d", files.size() );
        File tmp_file = File.createTempFile( prefix, suffix, tmp_root );
        
        tmp_file.deleteOnExit();
        files.add( tmp_file );
        return tmp_file;
    }
    
    
    public synchronized File
    spillMap( Map<Object, List<T1>> map )
    {
        long started = System.currentTimeMillis();
        try
        {
            File f = getTempFile();
            if( verbose )
            {
                System.out.printf( "spill: %s", f.getAbsolutePath() );
            }

            ObjectOutputStream oos = openOutputStream( f );
            //spill & count records
            int i = 0;
            for( Entry<Object, List<T1>> entry : map.entrySet() )
            {
                oos.writeObject( new Pair<Object, List<T1>>( entry.getKey(), entry.getValue() ) );
               
                for( T1 e : entry.getValue() )
                    if( null != e )
                        ++ i;
                
                entry.setValue( null );
                oos.reset();
            }

            oos.flush();
            oos.close();
           
            if( verbose )
            {
                System.out.printf( "\t-%d map entries ( %d records ) in %d ms\n", 
                                   map.size(), 
                                   i, 
                                   System.currentTimeMillis() - started );
            }
            return f;
        } catch( Exception e )
        {
            throw new RuntimeException( e );
        } 
    }


    private ObjectOutputStream 
    openOutputStream( File file ) throws IOException, FileNotFoundException
    {
        ObjectOutputStream oos = new ObjectOutputStream( new BufferedOutputStream( new GZIPOutputStream( new BufferedOutputStream( new FileOutputStream( file ) ), OUTPUT_BUFFER_SIZE ) 
                                                                                   { 
                                                                                       { 
                                                                                           def.setLevel( Deflater.BEST_SPEED ); 
                                                                                       } 
                                                                                   }, 
                                                                                   OUTPUT_BUFFER_SIZE  
                                                                                 ) 
                                                       );
        return oos;
    }

    
    public synchronized Map<Object, List<T1>>
    fillMap( File file )
    {
        Map<Object, List<T1>> result = new HashMap<Object, List<T1>>( spill_page_size );
        ObjectInputStream ois = null;
        long started = System.currentTimeMillis();
        int  i = 0;
        
        try
        {
            ois = openInputStream( file );
            if( verbose )
            {
                System.out.printf( "fill:  %s", file.getAbsolutePath() );
            }
            for( ; ; )
            {
                @SuppressWarnings( "unchecked" )
                Pair<Object, List<T1>> entry = (Pair<Object, List<T1>>)ois.readObject();
                result.put( entry.key, entry.value );
                for( T1 e : entry.value )
                    if( null != e )
                        ++ i;
                
            }
            
        } catch( EOFException eof )
        {
            if( null != ois )
                try
                {
                    ois.close();
                } catch( IOException e )
                {
                    e.printStackTrace();
                }
            if( verbose )
            {
                System.out.printf( "\t+%d map entries ( %d records ) in %d ms\n", 
                                   result.size(), 
                                   i, 
                                   System.currentTimeMillis() - started );
            }
            return result;
        } catch( Exception e )
        {
            throw new RuntimeException( e );
        } 
    }

    
    public List<T1>
    newListBucket() 
    {
        if( use_spill && spill_page_size <= super.objects.size() )
        {
            spillMap( super.objects );
            super.objects.clear();
        }
        
        return super.newListBucket();
    }
    
    
    public synchronized void
    cascadeErrors() throws DataEaterException
    {
        if( use_spill && files.size() > 0 )
        {
            try
            {
                use_spill = false;
                int generation = -1;
                int i = 0;
                do
                {
                    if( super.objects.isEmpty() )
                        super.objects = fillMap( files.get( i++ ) );    
                    
                    generation = files.size();
                    for( int j = i; j < generation;  ++ j )
                    {
                        ObjectInputStream ois = null;
                        ObjectOutputStream oos = null;
                        int read_spots = 0;
                        int save_spots = 0;
                        
                        
                        try
                        {
                            File f = files.get( j );
                            ois = openInputStream( f );
                            if( verbose )
                                System.out.printf( "fill:  %s", f.getAbsolutePath() );

                            for( ; ; )
                            {
                                @SuppressWarnings( "unchecked" )
                                Pair<Object, List<T1>> entry = (Pair<Object, List<T1>>)ois.readObject();
                                ++read_spots;
                                if( super.objects.containsKey( entry.key ) )
                                {
                                    for( T1 object : entry.value )
                                        if( null != object )
                                            eat( object );
                                } else
                                {
                                    if( null == oos )
                                        oos = openOutputStream( getTempFile() );

                                    oos.writeObject( entry );
                                    
                                    entry.value = null;
                                    oos.reset();
                                    
                                    ++save_spots;
                                }
                            }
                        } catch( EOFException eof )
                        {
                            if( verbose )
                            {
                                System.out.printf( "\t+%d records\n", 
                                                   read_spots );
                            }
                            ois.close();
                        }
                        
                        if( null != oos )
                        {
                            
                            if( verbose )
                            {
                                System.out.printf( "spill: %s\t-%d records\n", 
                                                   files.get( files.size() - 1 ).getAbsolutePath(), 
                                                   save_spots );
                            }
                            oos.close();
                        }
                    }
                    
                    // flush unassembled
                    super.cascadeErrors();
                    super.objects.clear();
                    
                } while( ( i = generation ) < files.size() );
            } catch( Exception e )
            {
                throw new RuntimeException( e ); 
            }
            
        } else
            super.cascadeErrors();
        
    }


    private 
    ObjectInputStream 
    openInputStream( File file ) throws IOException, FileNotFoundException
    {
        return new ObjectInputStream( new BufferedInputStream( new GZIPInputStream( new BufferedInputStream( new FileInputStream( file ) ) ) ) );
    }
}
