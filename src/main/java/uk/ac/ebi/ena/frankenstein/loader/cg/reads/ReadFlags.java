package uk.ac.ebi.ena.frankenstein.loader.cg.reads;

public enum 
ReadFlags
{
    _BothMatches,                // 0
    LeftHalfDnbNoMatches,        // 1
    LeftHalfDnbMapOverflow,      // 2

    RightHalfDnbNoMatches,       // 4
    _BothNoMatches,              // 5
    _LeftOverflowRightNoMatches, // 6

    RightHalfDnbMapOverflow,     // 8
    _LeftNoMatchesRightOverflow, // 9
    _BothOverflow;               //10
    
    
    public static ReadFlags
    forString( String number )
    {
        if( 1 == number.length() )
        {
            switch( number.charAt( 0 ) )
            {
            case '0': return _BothMatches;
            case '1': return LeftHalfDnbNoMatches;
            case '2': return LeftHalfDnbMapOverflow;                 
            
            case '4': return RightHalfDnbNoMatches;
            case '5': return _BothNoMatches;
            case '6': return _LeftOverflowRightNoMatches;
            case '8': return RightHalfDnbMapOverflow;
            case '9': return _LeftNoMatchesRightOverflow;

            default: 
                throw new IllegalArgumentException( number );
            }
            
            
        } else if( 2 == number.length() )
        {
            if( number.charAt( 0 ) == '1' 
                && number.charAt( 1 ) == '0' )
                return _BothOverflow;
            throw new IllegalArgumentException( number );
        } else
        {
            throw new IllegalArgumentException( number );
        }
    }
    
}
