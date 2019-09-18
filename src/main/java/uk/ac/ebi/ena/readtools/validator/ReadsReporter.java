package uk.ac.ebi.ena.readtools.validator;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.StandardOpenOption;
import java.util.List;

import uk.ac.ebi.ena.readtools.webin.cli.rawreads.ScannerMessage;
import uk.ac.ebi.ena.readtools.webin.cli.rawreads.ScannerMessage.ScannerErrorMessage;

public class 
ReadsReporter 
{

    public void 
    write( File file, Severity severity, String origin, String message )
    {
        origin  = ( null == origin )  ? "" : origin;
        message = ( null == message ) ? "" : message;
        
        try
        {
            String out = String.format( "%s: %s%s", severity.toString(), origin, message );
            Files.write( file.toPath(),
                         out.getBytes( StandardCharsets.UTF_8 ),
                         StandardOpenOption.APPEND, StandardOpenOption.SYNC, StandardOpenOption.CREATE );

        } catch( IOException e )
        {
            throw new RuntimeException( e );
        }
    }

    
    public void 
    write( File file, ScannerMessage scannerMessage )
    {
        Severity severity = scannerMessage instanceof ScannerErrorMessage ? Severity.ERROR : Severity.INFO;
        String msg = scannerMessage.getMessage();
        if( null != scannerMessage.getThrowable() )
        {
            msg += scannerMessage.getThrowable().toString();
        }
        write( file, severity, scannerMessage.getOrigin(), msg );
    }

    
    public void 
    write( File file, List<ScannerMessage> messages )
    {
        for( ScannerMessage sm : messages )
        {
            write( file, sm );
        }
    }

    
    public enum 
    Severity
    {
        FIX,
        INFO,
        WARNING,
        ERROR
    }
}