package uk.ac.ebi.ena.readtools.refactored.validator;

import uk.ac.ebi.ena.readtools.refactored.provider.ReadsProvider;
import uk.ac.ebi.ena.readtools.refactored.read.FastqRead;

import java.util.Iterator;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import static uk.ac.ebi.ena.readtools.refactored.validator.InsdcReadsValidator.IUPAC_CODES;

public class FastqReadsValidator implements ReadsValidator<FastqRead> {
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
    final static Pattern pCasava18Name = Pattern.compile(
            "^@(.+)( +|\\t+)([0-9]+):[YN]:[0-9]*[02468]($|:.*$)");

    // regexs
    final static Pattern pBaseName = Pattern.compile("^@(.*)"); // for name of the record
    final static private Pattern pBases = Pattern.compile("^([" + IUPAC_CODES + "]*?)\\+$"); // bases, trailing '+' is obligatory
    final static private Pattern pQuals = Pattern.compile("^([!-~]*?)$"); //qualities


    @Override
    public boolean validate(ReadsProvider<FastqRead> provider) throws ReadsValidationException {
        long readCount = 0;
        for (FastqRead read : provider) {
            readCount++;
            validateRead(read, readCount);
        }
        return true;
    }

    private void validateRead(FastqRead read, long readCount) throws ReadsValidationException {
        validateReadName(read.getName(), readCount);
        validateBases(read.getBases(), readCount);
        validateQualityScores(read.getQualityScores(), readCount);
    }

    private void validateReadName(String name, long readCount) throws ReadsValidationException {
        Matcher matcher = pCasava18Name.matcher(name);
        if (!matcher.matches()) {
            matcher = pBaseName.matcher(name);
            if (!matcher.matches()) {
                throw new ReadsValidationException("Invalid read name: " + name, readCount);
            }
        }
    }

    private void validateBases(String bases, long readCount) throws ReadsValidationException {
        Matcher matcher = pBases.matcher(bases);
        if (!matcher.matches()) {
            throw new ReadsValidationException("Invalid bases: " + bases, readCount);
        }
    }

    private void validateQualityScores(String qualityScores, long readCount) throws ReadsValidationException {
        Matcher matcher = pQuals.matcher(qualityScores);
        if (!matcher.matches()) {
            throw new ReadsValidationException("Invalid quality scores: " + qualityScores, readCount);
        }
    }
}
