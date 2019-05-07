package uk.ac.ebi.ena.readtools.loader.fastq;

import java.io.EOFException;
import java.io.IOException;
import java.io.InputStream;
import java.io.Serializable;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import uk.ac.ebi.ena.readtools.loader.common.QualityNormalizer;
import uk.ac.ebi.ena.readtools.loader.common.QualityNormalizer.QualityNormaizationException;
import uk.ac.ebi.ena.readtools.loader.common.feeder.DataFeederException;
import uk.ac.ebi.ena.readtools.loader.common.feeder.FeedableData;
import uk.ac.ebi.ena.readtools.loader.common.feeder.FeedableDataChecker;

public class 
DataSpot implements Serializable
{
    public static class 
    DataSpotParams
    {
        public Long      line_no;
        public boolean   allow_empty;
        public ReadStyle read_style;
        public Matcher   m_bases;
        public Matcher   m_qname;
        public Matcher   m_quals;
        public Matcher   m_quals_sd;
        public Matcher   m_casava_1_8_name;
        public Matcher   m_base_name;


        public 
        DataSpotParams( Long      line_no,
                        boolean   allow_empty,
                        ReadStyle read_style, 
                        Matcher   m_bases, 
                        Matcher   m_qname, 
                        Matcher   m_quals, 
                        Matcher   m_quals_sd, 
                        Matcher   m_casava_1_8_name, 
                        Matcher   m_base_name )
        {
            this.line_no      = line_no;
            this.allow_empty  = allow_empty;
            this.read_style   = read_style;
            this.m_bases      = m_bases;
            this.m_qname      = m_qname;
            this.m_quals      = m_quals;
            this.m_quals_sd   = m_quals_sd;
            this.m_casava_1_8_name = m_casava_1_8_name;
            this.m_base_name       = m_base_name;
        }
    }


    static final long serialVersionUID = 1L;
    static final protected char[] line_separator = System.getProperty ( "line.separator" ).toCharArray();
    
    @FeedableData( method = "feedBaseName" )
    public String bname; // name for bases 
        
    @FeedableData( method = "feedBases" )
    public String bases;// bases
        
    @FeedableData( method = "feedQualName" )
    public String qname; // name for qualities // no + here, it will be in stop symbols

    @FeedableData( method = "feedQuals" )
    public String quals;

    
    @FeedableDataChecker
    public void
    checkFeed() throws QualityNormaizationException, DataFeederException
    {
        if( !params.allow_empty )
        {
            if( null == bases || null == quals 
                || 0 == bases.length() || 0 == quals.length() )
                throw new DataFeederException( params.line_no, "Empty lines not allowed" );
        }
        normailzer.normalize( quals );
    }
    

    private String stream_key;
    transient private int expected_base_length;
    transient private int expected_qual_length;
    transient private QualityNormalizer normailzer;
    
    
    // default - no constraints
    public
    DataSpot()
    {
        expected_base_length = -1;
        expected_qual_length = -1;
        normailzer = QualityNormalizer.SANGER;
        stream_key = null;
    }
    
    
    // default - no constraints
    public
    DataSpot( QualityNormalizer normalizer, 
              String            default_stream_attr,
              DataSpotParams    params )
    {
        this();
        this.normailzer = normalizer;
        this.stream_key = default_stream_attr;
        this.params     = params;
    }
   
   
    
    public 
    DataSpot( int expected_length )
    {
        this();
        expected_base_length = expected_length;
        params     = defaultParams();
    }
    
    
    public String
    getStreamKey()
    {
        return stream_key;
    }
    
    
    /*
    @ Each sequence identifier line starts with @
1    <instrument> Characters
    allowed:
    a-z, A-Z, 0-9 and
    underscore
2    Instrument ID
    <run number> Numerical Run number on instrument
3    <flowcell
    ID>
    Characters
    allowed:
    a-z, A-Z, 0-9
4    <lane> Numerical Lane number
5    <tile> Numerical Tile number
6    <x_pos> Numerical X coordinate of cluster
7    <y_pos> Numerical Y coordinate of cluster
SPACE HERE
8    <read> Numerical Read number. 1 can be single read or read 2 of pairedend
9    <is
    filtered>
    Y or N Y if the read is filtered, N otherwise
10    <control
    number>
    Numerical 0 when none of the control bits are on, otherwise it is
    an even number. See below.
11    <index
    sequence>
    ACTG Index sequence
    */
    //                                                          1        :  2   :    3       :   4  :  5   :   6   :  7          8 :  9 :  10         : 11
//    final static private Pattern p_casava_1_8_name = Pattern.compile( "^@([a-zA-Z0-9_-]+:[0-9]+:[a-zA-Z0-9]+:[0-9]+:[0-9]+:[0-9-]+:[0-9-]+) ([12]):[YN]:[0-9]*[02468]:[ACGTN]+$" );
    // relexed regular expression 
    final static Pattern p_casava_1_8_name = Pattern.compile( "^@([a-zA-Z0-9_-]+:[0-9]+:[a-zA-Z0-9_-]+:[0-9]+:[0-9]+:[0-9-]+:[0-9-]+) ([12]):[YN]:[0-9]*[02468]:.*$" );
    
    // regexs
    final static private Pattern p_base_name = Pattern.compile( "^@(\\S*)( .*$|$)" ); // for name of the record
    final static private Pattern p_bases     = Pattern.compile( "^([ACGTNactgn.]*?)\\+$" ); // bases, trailing '+' is obligatory
    final static private Pattern p_qual_name = Pattern.compile( "^(\\S*)(?: .*$|$)" );  // name of quality record 
  //  Pattern p_quals     = Pattern.compile( "^([\\!\\\"\\#\\$\\%\\&\\'\\(\\)\\*\\+,\\-\\.\\/0-9:;<=>\\?\\@A-I]+)$" ); //qualities
    final static private Pattern p_quals_sd  = Pattern.compile( "^-?[0-9]{0,3}( +-?[0-9]{1,3})*[ -]*?$" );
    final static private Pattern p_quals     = Pattern.compile( "^([!-~]*?)$" ); //qualities
    final static private char base_stopper = '+';

    transient private DataSpotParams params = null;
          
                                                  
    public static DataSpotParams                                               
    defaultParams()
    {
        return new DataSpotParams( 1L, 
                                   true,
                                   (ReadStyle) null, 
                                   p_bases.matcher( "" ), 
                                   p_qual_name.matcher( "" ), 
                                   p_quals.matcher( "" ), 
                                   p_quals_sd.matcher( "" ), 
                                   p_casava_1_8_name.matcher( "" ), 
                                   p_base_name.matcher( "" ) );
    
    }
    

    public enum
    ReadStyle
    {
        FASTQ,
        CASAVA18
    }
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //feeding methods. NOTE: always use "public" access modifier because of method accession performance issues!
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    public void 
    feedBaseName( InputStream is ) throws IOException, DataFeederException
    {
        String line = readLine( is );
        while( line.trim().length() == 0 )
            line = readLine( is );
        
        //Try to determine read style
        if( null == params.read_style )
        {
            //try casava 1.8
            if( params.m_casava_1_8_name.reset( line ).find() )
                params.read_style = ReadStyle.CASAVA18; 
            else
                params.read_style = ReadStyle.FASTQ;
        }
        
        switch( params.read_style )
        {
            case CASAVA18:
                if( !params.m_casava_1_8_name.reset( line ).find() )
                    throw new DataFeederException( params.line_no, String.format( "Line [%s] does not match %s regexp", line, ReadStyle.CASAVA18 ) );
                bname = String.format( "%s/%s", params.m_casava_1_8_name.group( 1 ), params.m_casava_1_8_name.group( 2 ) );
                break;
            
            case FASTQ:
                if( !params.m_base_name.reset( line ).find() )
                    throw new DataFeederException( params.line_no, String.format( "Line [%s] does not match %s regexp", line, ReadStyle.FASTQ ) );
                
                bname = params.m_base_name.group( 1 );
                break;
                
            default:
                throw new DataFeederException( params.line_no, String.format( "Line [%s] has no read style defined!", line ) );
        }
    }
    
    
    // get bases
    public void 
    feedBases( InputStream is ) throws IOException, DataFeederException
    {
        String line = readLine( is, -1, base_stopper );
        while( line.trim().length() == 0 )
            line = readLine( is, -1, base_stopper );

        if( !params.m_bases.reset( line ).find() )
            throw new DataFeederException( params.line_no, String.format( "Line [%s] does not match regexp", line ) );
            
        String value = params.m_bases.group( 1 );

        //check against expected, if any
        if( -1 < expected_base_length 
            && expected_base_length != value.length() )
            throw new DataFeederException( params.line_no, String.format( "Expected base length [%d] does not match readed one[%d]", expected_base_length, value.length() ) );    
                
        expected_qual_length = 0 == value.length() ? -1 : value.length();
        bases = value;
    }
    
    
    // get name of quality line
    public void 
    feedQualName( InputStream is ) throws IOException, DataFeederException
    {
        String line = readLine( is );        
        if( !params.m_qname.reset( line ).find() )
            throw new DataFeederException( params.line_no, String.format( "Line [%s] does not match regexp", line ) );
        qname = params.m_qname.group( 1 );
    }

    
    public void 
    feedQuals( InputStream is ) throws IOException, QualityNormaizationException, DataFeederException
    {
        String line = readLine( is, expected_qual_length );
        while( expected_qual_length >= 0 && line.trim().length() == 0 )
            line = readLine( is, expected_qual_length );

        String value = null;
        
        if( !params.m_quals.reset( line ).find() )
        {
            if( !params.m_quals_sd.reset( line ).matches() )
                throw new DataFeederException( params.line_no, String.format( "Line [%s] does not match regexp", line ) );
            else
            {
                line += readLine( is );

                if( !params.m_quals_sd.reset( line ).matches() )
                    throw new DataFeederException( params.line_no, String.format( "Line [%s] does not match regexp", line ) );
                
                String[] scores = line.split( " +" );
                StringBuilder sb = new StringBuilder( scores.length );
                for( String score : scores )
                    sb.append( (char) ( Integer.parseInt( score ) + '!' ) );
                value = sb.toString();
                
                if( expected_qual_length != value.length() )
                    throw new DataFeederException( params.line_no, String.format( "%s Expected qual length [%d] does not match length of readed one[%d]", 
                                                                                 bname, 
                                                                                 expected_qual_length, 
                                                                                 value.length() ) );    

            }
        } else
        {
            value = params.m_quals.group( 1 );
        
            //check against expected
            if( expected_qual_length >= 0 && expected_qual_length != value.length() )
                throw new DataFeederException( params.line_no, String.format( "%s Expected qual length [%d] does not match length of readed one[%d]", 
                                                                             bname, 
                                                                             expected_qual_length, 
                                                                             value.length() ) );    

            //we are lenient now.
            if( expected_qual_length >= 0 )
            {
                try
                {
                    line = readLine( is );
                    if( line.trim().length() > 0 )
                        throw new DataFeederException( params.line_no, String.format( "Trailing character(s) [%s] after expected number of quals [%d]", 
                                                                                     line, 
                                                                                    expected_qual_length ) );
                } catch( IOException e )
                {
                    ;
                }
            }
        }
        quals = value;
    }

    
    // reads stream line till line separator
    protected String
    readLine( InputStream istream ) throws IOException
    {
        return readLine( istream, -1 );
    }

    
    // reads stream line till either the length or in case of len = -1 till line separator 
    protected String
    readLine( InputStream istream, 
              long len ) throws IOException
    {
            return readLine( istream, len, -1 );
    }
        

    // reads stream line till either the length or in case of len = -1 till line separator
    // here is unfair checking of stop_seq ( only 1 symbol is allowed!!! )
    // TODO - split and simplyfy
    protected String
    __readLine( InputStream istream, 
                int       len,
                int       stop ) throws IOException
    {
        assert istream.markSupported();
        StringBuilder b = new StringBuilder( 1024 );
        
        int space_cnt = 0;
out:    for( int i = 0; len == -1 || i < len; )
        {
            //try to read at least len bytes
            byte buf[] = new byte[ 1024 ];
            int read_len = -1 == len ? buf.length : len - i;
            istream.mark( read_len );
            read_len = istream.read( buf, 0, read_len );
            if( -1 == read_len )
                throw new EOFException();
            for( int real_i = 0; real_i < read_len; ++ i, ++real_i )
            {
                int c = buf[ real_i ];
                if( c == ' ' )
                {
                    space_cnt++;
                } else if( c == '\r' )
                {
                    --i;
                    continue;
                } else if( c == '\n' )
                {   
                    params.line_no++;
                    //get rid of trailing spaces;
                    if( space_cnt > 0 )
                    {
                        b.delete( b.length() - space_cnt, b.length() );
                        i -= space_cnt;
                        space_cnt -= space_cnt;
                    }
                    
                    //return if no actual stop symbol and len is infinite
                    if( -1 == stop && -1 == len ) 
                    {
                        istream.reset();
                        istream.skip( real_i + 1 );
                        break out;
                    }
                    --i;
                    continue;
                } else if( c == -1 )
                {
                    //TODO: should remove trailing spaces?
                    if( b.length() > 0 )
                    {
                        istream.reset();
                        istream.skip( real_i + 1 );
                        return b.toString();
                    }
                } else
                {
                    space_cnt -= space_cnt;
                }
                
                b.append( (char)c );
                if( c == stop )
                {
                    istream.reset();
                    istream.skip( real_i + 1 );
                    break out;
                }
            }
        } 
        return b.toString();
    }

    
    // reads stream line till either the length or in case of len = -1 till line separator
    // here is unfair checking of stop_seq ( only 1 symbol is allowed!!! )
    protected String
    readLine( InputStream istream, 
              long          len,
              int           stop ) throws IOException
    {
        StringBuilder b = new StringBuilder( 1024 );
        
        int space_cnt = 0;
        
        for( int i = 0; len == -1 || i < len; ++ i )
        {
            int c = istream.read();
            if( c == ' ' )
            {
                space_cnt++;
            } else if( c == '\r' )
            {
                --i;
                continue;
            } else if( c == '\n' )
            {   
                params.line_no++;
                //get rid of trailing spaces;
                if( space_cnt > 0 )
                {
                    b.delete( b.length() - space_cnt, b.length() );
                    i -= space_cnt;
                    space_cnt -= space_cnt;
                }
                
                //return if no actual stop symbol and len is infinite
                if( -1 == stop && -1 == len ) 
                    break;
                --i;
                continue;
            } else if( c == -1 )
            {
                //TODO: should remove trailing spaces?
                if( b.length() > 0 )
                    return b.toString();
                
                throw new EOFException();        
            } else
            {
                space_cnt -= space_cnt;
            }
            
            b.append( (char)c );
            if( c == stop )
                break;
        } 
        return b.toString();
    }

    
    public String
    toString()
    {
        return new StringBuilder()
            .append( "base_name = [" )
            .append( bname )
            .append( "]\n" )
            .append( "bases = [" )
            .append( bases )
            .append( "], length = " )
            .append( bases.length() )
            .append( "\nqual_name = [" )
            .append( qname )
            .append( "]\n" )
            .append( "quals = [" )
            .append( quals )
            .append( "], length = " )
            .append( quals.length() )
            .toString();
    }
}
