package uk.ac.ebi.ena.readtools.validator;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.stream.Collectors;

import htsjdk.samtools.SAMFormatException;
import htsjdk.samtools.cram.CRAMException;
import uk.ac.ebi.ena.readtools.validator.ReadsReporter.Severity;
import uk.ac.ebi.ena.readtools.webin.cli.rawreads.BamScanner;
import uk.ac.ebi.ena.readtools.webin.cli.rawreads.FastqScanner;
import uk.ac.ebi.ena.readtools.webin.cli.rawreads.RawReadsFile;
import uk.ac.ebi.ena.readtools.webin.cli.rawreads.RawReadsFile.Filetype;
import uk.ac.ebi.ena.readtools.webin.cli.rawreads.ScannerMessage;
import uk.ac.ebi.ena.readtools.webin.cli.rawreads.ScannerMessage.ScannerErrorMessage;
import uk.ac.ebi.ena.readtools.webin.cli.rawreads.refs.CramReferenceInfo;
import uk.ac.ebi.ena.webin.cli.validator.api.ValidationResponse.status;
import uk.ac.ebi.ena.webin.cli.validator.api.Validator;
import uk.ac.ebi.ena.webin.cli.validator.file.SubmissionFile;
import uk.ac.ebi.ena.webin.cli.validator.manifest.ReadsManifest;
import uk.ac.ebi.ena.webin.cli.validator.manifest.ReadsManifest.FileType;
import uk.ac.ebi.ena.webin.cli.validator.response.ReadsValidationResponse;

public class ReadsValidator 
implements Validator<ReadsManifest, ReadsValidationResponse> 
{

    private ReadsManifest manifest;
    private ReadsReporter reporter = new ReadsReporter();

    @Override public ReadsValidationResponse 
    validate( ReadsManifest manifest )
    {
        this.manifest = manifest;
        if( manifest == null )
        {
            throw new RuntimeException( "Manifest can not be null." );
        }
        if( manifest.getReportFile() == null )
        {
            throw new RuntimeException( "Report file is missing." );
        }
        if( manifest.getProcessDir() == null )
        {
            reporter.write( manifest.getReportFile(), Severity.ERROR, "", "Process directory is missing." );
            throw new RuntimeException( "Process directory is missing." );
        }

        return validateSubmissionForContext();
    }

    
    private ReadsValidationResponse 
    validateSubmissionForContext()
    {
        boolean valid = true;
        AtomicBoolean paired = new AtomicBoolean();

        List<RawReadsFile> files = createReadFiles();

        for( RawReadsFile rf : files )
        {
            if( Filetype.fastq.equals( rf.getFiletype() ) )
            {
                valid = readFastqFile( files, paired );
            } else if( Filetype.bam.equals( rf.getFiletype() ) )
            {
                valid = readBamFile( files, paired );
            } else if( Filetype.cram.equals( rf.getFiletype() ) )
            {
                valid = readCramFile( files, paired );
            } else
            {
                throw new RuntimeException( "Unsupported file type: " + rf.getFiletype().name() );
            }
            break;
        }

        ReadsValidationResponse resp = new ReadsValidationResponse();
        resp.setStatus( valid ? status.VALIDATION_SUCCESS : status.VALIDATION_ERROR );
        resp.setPaired( paired.get() );
        return resp;
    }

    
    private boolean 
    readCramFile( List<RawReadsFile> files, AtomicBoolean paired )
    {
        boolean valid = true;
        CramReferenceInfo cri = new CramReferenceInfo();
        for( RawReadsFile rf : files )
        {
            try
            {
                Map<String, Boolean> ref_set = cri.confirmFileReferences( new File( rf.getFilename() ) );
                if( !ref_set.isEmpty() && ref_set.containsValue( Boolean.FALSE ) )
                {
                    reporter.write( rf.getReportFile().toFile(), Severity.ERROR,
                            "",
                            "Unable to find reference sequence(s) from the CRAM reference registry: " +
                                    ref_set.entrySet()
                                            .stream()
                                            .filter( e -> !e.getValue() )
                                            .map( e -> e.getKey() )
                                            .collect( Collectors.toList() ) );
                    valid = false;
                }

            } catch( IOException ioe )
            {
                reporter.write( rf.getReportFile().toFile(), Severity.ERROR,
                        "",
                        ioe.getMessage() ); // RuntimeEx
                valid = false;
            }
        }

        return valid && readBamFile( files, paired );
    }

    
    private boolean 
    readFastqFile( List<RawReadsFile> files, AtomicBoolean paired )
    {
        try
        {
            if( files.size() > 2 )
            {
                String msg = "Unable to validate unusual amount of files: " + files;
                reportToFileList( files, msg );
                return false;
            }

            FastqScanner fs = new FastqScanner( manifest.getPairingHorizon() ) 
            {
                @Override protected void 
                logProcessedReadNumber( long count )
                {
                    ReadsValidator.this.logProcessedReadNumber( count );
                }
                
                
                @Override protected void 
                logFlushMsg( String msg )
                {
                    ReadsValidator.this.logFlushMsg( msg );
                    
                }
            };

            List<ScannerMessage> list = fs.checkFiles( files.toArray( new RawReadsFile[files.size()] ) );

            paired.set( fs.getPaired() );
            files.forEach( rf -> {
                reporter.write( rf.getReportFile().toFile(), list );
            } );

            return isValid( list );
        } catch( Exception ex )
        {
            throw new RuntimeException( ex );
        } catch( Throwable throwable )
        {
            throw new RuntimeException( throwable );
        }
    }

    
    private void 
    reportToFileList( List<RawReadsFile> files, String msg )
    {
        for( RawReadsFile rf : files )
        {
            reporter.write( rf.getReportFile().toFile(), Severity.ERROR, "", msg );
        }
    }

    
    private boolean 
    readBamFile( List<RawReadsFile> files, AtomicBoolean paired )
    {
        BamScanner scanner = new BamScanner()
        {
            @Override protected void 
            logProcessedReadNumber( long count )
            {
                ReadsValidator.this.logProcessedReadNumber( count );
            }
        };


        boolean valid = true;
        for( RawReadsFile rf : files )
        {
            try
            {
                String msg = String.format( "Processing file %s\n", rf.getFilename() );
                reporter.write( manifest.getReportFile(), Severity.INFO, "", msg );

                List<ScannerMessage> list = Filetype.cram == rf.getFiletype() ? scanner.readCramFile( rf, paired )
                        : scanner.readBamFile( rf, paired );
                list.stream().forEachOrdered( m -> reporter.write( rf.getReportFile().toFile(), m ) );
                valid = isValid( list );

            } catch( SAMFormatException | CRAMException e )
            {
                reporter.write( rf.getReportFile().toFile(), Severity.ERROR, "", e.getMessage() );
                valid = false;

            } catch( IOException ex )
            {
                throw new RuntimeException( ex );
            }
        }
        return valid;
    }


    private List<RawReadsFile> 
    createReadFiles()
    {
        ReadsManifest manifest = this.manifest;

        RawReadsFile.AsciiOffset asciiOffset = null;
        RawReadsFile.QualityScoringSystem qualityScoringSystem = null;

        if( manifest.getQualityScore() != null )
        {
            switch( manifest.getQualityScore() )
            {
            case PHRED_33:
                asciiOffset = RawReadsFile.AsciiOffset.FROM33;
                qualityScoringSystem = RawReadsFile.QualityScoringSystem.phred;
                break;
            case PHRED_64:
                asciiOffset = RawReadsFile.AsciiOffset.FROM64;
                qualityScoringSystem = RawReadsFile.QualityScoringSystem.phred;
                break;
            case LOGODDS:
                asciiOffset = null;
                qualityScoringSystem = RawReadsFile.QualityScoringSystem.log_odds;
                break;
            }
        }
        List<RawReadsFile> files = this.manifest.files().get()
                .stream()
                .map( file -> createRawReadsFile( file ) )
                .collect( Collectors.toList() );

        // Set FASTQ quality scoring system and ascii offset.

        for( RawReadsFile f : files )
        {
            if( f.getFiletype().equals( Filetype.fastq ) )
            {
                if( qualityScoringSystem != null )
                {
                    f.setQualityScoringSystem( qualityScoringSystem );
                }
                if( asciiOffset != null )
                {
                    f.setAsciiOffset( asciiOffset );
                }
            }
        }

        return files;
    }

    
    public static RawReadsFile 
    createRawReadsFile( SubmissionFile<FileType> file )
    {
        Path inputDir = file.getFile().toPath().getParent();

        RawReadsFile f = new RawReadsFile();
        f.setInputDir( inputDir );
        f.setReportFile( file.getReportFile().toPath() );
        f.setFiletype( Filetype.valueOf( file.getFileType().name().toLowerCase() ) );

        String fileName = file.getFile().getPath();

        if( !Paths.get( fileName ).isAbsolute() )
        {
            f.setFilename( inputDir.resolve( Paths.get( fileName ) ).toString() );
        } else
        {
            f.setFilename( fileName );
        }

        return f;
    }

    
    private void 
    logProcessedReadNumber( long count )
    {
        String msg = String.format( "\rProcessed %16d read(s)", count );
        logFlushMsg( msg );
    }

    
    private void
    logFlushMsg( String msg )
    {
        System.out.print( msg );
        System.out.flush();
    }
    
    
    boolean 
    isValid( List<ScannerMessage> list )
    {
        return !list.stream().anyMatch( e -> e instanceof ScannerErrorMessage );
    }
}
