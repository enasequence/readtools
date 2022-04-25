package uk.ac.ebi.ena.readtools.fastq.ena;

import org.apache.commons.exec.CommandLine;
import org.apache.commons.exec.DefaultExecutor;
import org.apache.commons.exec.Executor;
import org.apache.commons.exec.PumpStreamHandler;

import java.io.IOException;
import java.io.InputStream;
import java.io.PipedInputStream;
import java.io.PipedOutputStream;
import java.nio.file.Path;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

public class Sra2SamInputStream extends InputStream {
    private PipedInputStream sra2samStream;
    private Future<Integer> future;
    private Path filePath;

    public Sra2SamInputStream(Path filePath) throws IOException {
        this.filePath = filePath;

        String cmd = "bam-dump";
        CommandLine commandLine = new CommandLine(cmd);
        commandLine.addArgument(filePath.toAbsolutePath().toString(), false);

        PipedOutputStream stdoutStream = new PipedOutputStream();

        sra2samStream = new PipedInputStream();
        stdoutStream.connect(sra2samStream);

        Executor apacheExecutor = new DefaultExecutor();
        apacheExecutor.setExitValues(null);
        apacheExecutor.setStreamHandler(new PumpStreamHandler(stdoutStream, System.err));

        ExecutorService executor = Executors.newSingleThreadExecutor();

        future = executor.submit(() -> apacheExecutor.execute(commandLine));
    }

    @Override
    public int read() throws IOException {
        return sra2samStream.read();
    }

    @Override
    public void close() throws IOException {
        sra2samStream.close();
        try {
            if (!future.isDone()) {
                future.cancel(true);
            }
            Integer exitCode = future.get();
            if (!exitCode.equals(0)) {
                throw new IOException("Failed to read " + filePath + " bam-dump exit code " + exitCode);
            }
        } catch (ExecutionException | InterruptedException e) {
            throw new IOException(e);
        }
    }
}
