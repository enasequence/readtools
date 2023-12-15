package uk.ac.ebi.ena.readtools.refactored;

import org.junit.Assert;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import java.util.zip.GZIPOutputStream;

public class TestFileUtil {
    public static File
    createOutputFolder() throws IOException {
        File output = File.createTempFile( "test", "test" );
        Assert.assertTrue( output.delete() );
        Assert.assertTrue( output.mkdirs() );
        return output;
    }

    public static Path
    saveRandomized(String content, Path folder, boolean gzip, String... suffix) throws IOException {
        Path file = Files.createTempFile( "_content_", "_content_" );
        Files.write( file, content.getBytes( StandardCharsets.UTF_8 ), StandardOpenOption.TRUNCATE_EXISTING, StandardOpenOption.SYNC );
        Path path = Files.createTempFile( folder, "COPY", ( file.getName( file.getNameCount() - 1 ) + ( suffix.length > 0 ? Stream.of( suffix ).collect( Collectors.joining( ".", ".", "" ) ) : "" ) ));
        OutputStream os;
        Files.copy( file, ( os = gzip ? new GZIPOutputStream( Files.newOutputStream( path, StandardOpenOption.TRUNCATE_EXISTING, StandardOpenOption.SYNC ) )
                : Files.newOutputStream( path, StandardOpenOption.TRUNCATE_EXISTING, StandardOpenOption.SYNC ) ) );
        os.flush();
        os.close();
        Assert.assertTrue( Files.exists( path ) );
        Assert.assertTrue( Files.isRegularFile( path ) );
        return path;
    }
}
