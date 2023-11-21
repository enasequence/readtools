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
package uk.ac.ebi.ena.readtools.loader.common.converter;

import java.io.EOFException;
import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.util.Set;
import java.util.function.Consumer;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import uk.ac.ebi.ena.readtools.common.reads.QualityNormalizer;
import uk.ac.ebi.ena.readtools.loader.common.InvalidBaseCharacterException;
import uk.ac.ebi.ena.readtools.loader.fastq.Read;

class ReadReader {

    public enum ReadStyle {
        FASTQ,
        CASAVA18
    }

    public static class DataSpotReaderParams
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


        public DataSpotReaderParams(Long      line_no,
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
    // relaxed regular expression
    final static Pattern p_casava_1_8_name = Pattern.compile(
            "^@(.+)( +|\\t+)([0-9]+):[YN]:[0-9]*[02468]($|:.*$)" );

    // regexs
    final static Pattern p_base_name = Pattern.compile( "^@(.*)" ); // for name of the record
    final static private Pattern p_bases     = Pattern.compile( "^([ACGTNUactgnu.]*?)\\+$" ); // bases, trailing '+' is obligatory
    final static private Pattern p_qual_name = Pattern.compile( "^(\\S*)(?: .*$|$)" );  // name of quality record
    //  Pattern p_quals     = Pattern.compile( "^([\\!\\\"\\#\\$\\%\\&\\'\\(\\)\\*\\+,\\-\\.\\/0-9:;<=>\\?\\@A-I]+)$" ); //qualities
    final static private Pattern p_quals_sd  = Pattern.compile( "^-?[0-9]{0,3}( +-?[0-9]{1,3})*[ -]*?$" );
    final static private Pattern p_quals     = Pattern.compile( "^([!-~]*?)$" ); //qualities
    final static private char base_stopper = '+';

    private final DataSpotReaderParams params = defaultParams();

    /**
     * Default read index assigned to each {@link Read} spot returned by this reader instance.
     */
    private final String defaultReadIndex;

    private final Consumer<Read> readQualityNormalizer;

    private final int expectedBaseLength;

    /**
     * Creates a reader that does not perform any quality normalization.
     */
    public ReadReader(String defaultReadIndex) {
        this(null, defaultReadIndex);
    }

    public ReadReader(QualityNormalizer normalizer) {
        this(normalizer, null);
    }

    public ReadReader(QualityNormalizer normalizer, String defaultReadIndex ) {
        this(normalizer, defaultReadIndex, -1);
    }

    public ReadReader(QualityNormalizer normalizer, int expectedLength ) {
        this(normalizer, null, expectedLength);
    }

    public ReadReader(QualityNormalizer normalizer, String defaultReadIndex, int expectedLength ) {
        if (normalizer == null) {
            readQualityNormalizer = read -> {}; //No normalization is performed.
        } else {
            readQualityNormalizer = read -> {
                byte[] quals = read.getQualityScores().getBytes(StandardCharsets.UTF_8);
                normalizer.normalize(quals);
                read.quals = new String(quals, StandardCharsets.UTF_8);
            };
        }

        this.defaultReadIndex = defaultReadIndex;
        this.expectedBaseLength = expectedLength;
    }

    public Read read(InputStream inputStream) throws IOException {
        boolean recordStarted = false;

        Read read = new Read();
        read.defaultReadIndex = defaultReadIndex;

        try {
            readBaseName(inputStream, read);

            recordStarted = true;

            readBases(inputStream, read);
            readQualName(inputStream, read);
            readQuals(inputStream, read);

            checkReadData(read);

            return read;
        } catch (EOFException e) {
            // Read data must always be checked even when the end of stream has been reached.
            if( recordStarted ) {
                checkReadData(read);
            }

            throw e;
        }
    }

    private void readBaseName(InputStream is, Read read) throws IOException
    {
        String line = readLine( is );
        while( line.trim().length() == 0 )
            line = readLine( is );

        if (line.trim().length() > 257) // uint8 + @ character
            throw new ConverterException(
                    params.line_no,
                    String.format("Line's length exceeds 256 characters: [%s]", line));

        //Try to determine read style
        if( null == params.read_style )
        {
            //try casava 1.8
            if( params.m_casava_1_8_name.reset( line ).find() )
                params.read_style = ReadReader.ReadStyle.CASAVA18;
            else
                params.read_style = ReadReader.ReadStyle.FASTQ;
        }

        switch( params.read_style )
        {
            case CASAVA18:
                if( !params.m_casava_1_8_name.reset( line ).find() )
                    throw new ConverterException( params.line_no, String.format( "Line [%s] does not match %s regexp", line, ReadReader.ReadStyle.CASAVA18 ) );
                read.name = params.m_casava_1_8_name.group( 1 ) + "/" + params.m_casava_1_8_name.group( 3 );
                break;

            case FASTQ:
                if( !params.m_base_name.reset( line ).find() )
                    throw new ConverterException( params.line_no, String.format( "Line [%s] does not match %s regexp", line, ReadReader.ReadStyle.FASTQ ) );

                read.name = params.m_base_name.group( 1 );
                break;

            default:
                throw new ConverterException( params.line_no, String.format( "Line [%s] has no read style defined!", line ) );
        }
    }

    // get bases
    private void readBases(InputStream is, Read read) throws IOException
    {
        String line = readLine( is, -1, base_stopper );
        while( line.trim().length() == 0 )
            line = readLine( is, -1, base_stopper );

        if( !params.m_bases.reset( line ).find() )
            handleInvalidBases(line);

        String value = params.m_bases.group( 1 );

        //check against expected, if any
        if( -1 < expectedBaseLength
                && expectedBaseLength != value.length() )
            throw new ConverterException( params.line_no, String.format( "Expected base length [%d] does not match the read one[%d]", expectedBaseLength, value.length() ) );

        read.bases = value;
    }

    // get name of quality line
    private void readQualName(InputStream is, Read read) throws IOException
    {
        String line = readLine( is );
        if( !params.m_qname.reset( line ).find() )
            throw new ConverterException( params.line_no, String.format( "Line [%s] does not match regexp", line ) );
    }

    private void readQuals(InputStream is, Read read) throws IOException
    {
        int expectedQualLength = 0 == read.getBases().length() ? -1 : read.getBases().length();

        String line = readLine( is, expectedQualLength );
        while( expectedQualLength >= 0 && line.trim().length() == 0 )
            line = readLine( is, expectedQualLength );

        String value = null;

        if( !params.m_quals.reset( line ).find() )
        {
            if( !params.m_quals_sd.reset( line ).matches() )
                throw new ConverterException( params.line_no, String.format( "Line [%s] does not match regexp", line ) );
            else
            {
                line += readLine( is );

                if( !params.m_quals_sd.reset( line ).matches() )
                    throw new ConverterException( params.line_no, String.format( "Line [%s] does not match regexp", line ) );

                String[] scores = line.split( " +" );
                StringBuilder sb = new StringBuilder( scores.length );
                for( String score : scores )
                    sb.append( (char) ( Integer.parseInt( score ) + '!' ) );
                value = sb.toString();

                if( expectedQualLength != value.length() )
                    throw new ConverterException( params.line_no, String.format( "%s Expected qual length [%d] does not match length of the read one[%d]",
                            read.getName(),
                            expectedQualLength,
                            value.length() ) );

            }
        } else
        {
            value = params.m_quals.group( 1 );

            //check against expected
            if( expectedQualLength >= 0 && expectedQualLength != value.length() )
                throw new ConverterException( params.line_no, String.format( "%s Expected qual length [%d] does not match length of the read one[%d]",
                        read.getName(),
                        expectedQualLength,
                        value.length() ) );

            //we are lenient now.
            if( expectedQualLength >= 0 )
            {
                try
                {
                    line = readLine( is );
                    if( line.trim().length() > 0 )
                        throw new ConverterException( params.line_no, String.format(
                            "Found extra character(s) after expected number (%d) of quals. Extra characters : [%s]",
                            expectedQualLength, line) );
                } catch( IOException e )
                {
                    ;
                }
            }
        }
        read.quals = value;
    }

    private void checkReadData(Read read)
    {
        if( !params.allow_empty )
        {
            if( null == read.getBases() || null == read.getQualityScores()
                    || 0 == read.getBases().length() || 0 == read.getQualityScores().length() )
                throw new ConverterException( params.line_no, "Empty lines not allowed" );
        }

        readQualityNormalizer.accept(read);
    }

    // reads stream line till line separator
    protected String readLine( InputStream istream ) throws IOException
    {
        return readLine( istream, -1 );
    }

    // reads stream line till either the length or in case of len = -1 till line separator
    protected String readLine( InputStream istream, long len ) throws IOException
    {
        return readLine( istream, len, -1 );
    }

    // reads stream line till either the length or in case of len = -1 till line separator
    // here is unfair checking of stop_seq ( only 1 symbol is allowed!!! )
    protected String readLine( InputStream istream, long len, int stop ) throws IOException
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

    private void handleInvalidBases(String bases) {
        Matcher matcher = p_bases.matcher("");

        //remove the trailing '+' symbol.
        bases = bases.substring(0, bases.length() - 1);

        Set<Character> invalidBaseChars = bases.chars()
                .mapToObj(intChar -> (char)intChar)
                // '+' symbol is added because its present in the regex pattern.
                .filter(base -> !matcher.reset(String.valueOf(base) + "+").matches())
                .collect(Collectors.toSet());

        throw new InvalidBaseCharacterException(bases, invalidBaseChars);
    }

    private static DataSpotReaderParams defaultParams() {
        return new DataSpotReaderParams( 0L,
                true,
                (ReadStyle) null,
                p_bases.matcher( "" ),
                p_qual_name.matcher( "" ),
                p_quals.matcher( "" ),
                p_quals_sd.matcher( "" ),
                p_casava_1_8_name.matcher( "" ),
                p_base_name.matcher( "" ) );

    }
}
