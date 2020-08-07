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
import uk.ac.ebi.ena.readtools.old.webin.cli.rawreads.BamScanner;
import uk.ac.ebi.ena.readtools.old.webin.cli.rawreads.FastqScanner;
import uk.ac.ebi.ena.readtools.old.webin.cli.rawreads.RawReadsFile;
import uk.ac.ebi.ena.readtools.old.webin.cli.rawreads.RawReadsFile.Filetype;
import uk.ac.ebi.ena.readtools.old.webin.cli.rawreads.refs.CramReferenceInfo;
import uk.ac.ebi.ena.webin.cli.validator.api.ValidationResponse.status;
import uk.ac.ebi.ena.webin.cli.validator.api.Validator;
import uk.ac.ebi.ena.webin.cli.validator.file.SubmissionFile;
import uk.ac.ebi.ena.webin.cli.validator.manifest.ReadsManifest;
import uk.ac.ebi.ena.webin.cli.validator.manifest.ReadsManifest.FileType;
import uk.ac.ebi.ena.webin.cli.validator.message.ValidationMessage;
import uk.ac.ebi.ena.webin.cli.validator.message.ValidationOrigin;
import uk.ac.ebi.ena.webin.cli.validator.message.ValidationResult;
import uk.ac.ebi.ena.webin.cli.validator.response.ReadsValidationResponse;

public class ReadsValidator 
implements Validator<ReadsManifest, ReadsValidationResponse> 
{

    private ReadsManifest manifest;

    @Override public ReadsValidationResponse 
    validate( ReadsManifest manifest )
    {
        this.manifest = manifest;
        if( manifest == null )
        {
            throw new RuntimeException( "Manifest is missing." );
        }
        if( manifest.getReportFile() == null )
        {
            throw new RuntimeException( "Report file is missing." );
        }
        if( manifest.getProcessDir() == null )
        {
            throw new RuntimeException( "Process directory is missing." );
        }
        return validateSubmissionForContext();
    }

    private ReadsValidationResponse 
    validateSubmissionForContext()
    {
        ValidationResult result = new ValidationResult(manifest.getReportFile());

        AtomicBoolean paired = new AtomicBoolean();

        List<RawReadsFile> files = createReadFiles();

        // TODO: remove for loop
        for( RawReadsFile rf : files )
        {
            if( Filetype.fastq.equals( rf.getFiletype() ) )
            {
                readFastqFile( result, files, paired );
            } else if( Filetype.bam.equals( rf.getFiletype() ) )
            {
                readBamFile( result, files, paired );
            } else if( Filetype.cram.equals( rf.getFiletype() ) )
            {
                 readCramFile( result, files, paired );
            } else
            {
                throw new RuntimeException( "Unsupported file type: " + rf.getFiletype().name() );
            }
            break;
        }

        ReadsValidationResponse resp = new ReadsValidationResponse();
        resp.setStatus( result.isValid() ? status.VALIDATION_SUCCESS : status.VALIDATION_ERROR );
        resp.setPaired( paired.get() );
        return resp;
    }

    
    private void
    readCramFile( ValidationResult result, List<RawReadsFile> files, AtomicBoolean paired )
    {
        CramReferenceInfo cri = new CramReferenceInfo();
        for( RawReadsFile rf : files )
        {
            ValidationResult fileResult = result.create(new ValidationOrigin("file", rf.getFilename()));

            try
            {
                Map<String, Boolean> ref_set = cri.confirmFileReferences( new File( rf.getFilename() ) );
                if( !ref_set.isEmpty() && ref_set.containsValue( Boolean.FALSE ) )
                {
                    fileResult.add(ValidationMessage.error(
                                    "Unable to find reference sequence(s) from the CRAM reference registry: " +
                                    ref_set.entrySet()
                                            .stream()
                                            .filter( e -> !e.getValue() )
                                            .map( e -> e.getKey() )
                                            .collect( Collectors.toList() ) ));
                }

            } catch( IOException ex )
            {
                fileResult.add(ValidationMessage.error(ex));
            }
        }

        if (!result.isValid()) {
            return;
        }

        readBamOrCramFile(result, files, paired);
    }
    
    private void
    readFastqFile( ValidationResult result, List<RawReadsFile> files, AtomicBoolean paired )
    {
        try
        {
            if( files.size() > 2 )
            {
                result.add(ValidationMessage.error("More than two fastq files were provided: " + files.size()));
                return;
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

            fs.checkFiles( result, files.toArray( new RawReadsFile[files.size()] ));
            paired.set( fs.getPaired() );
        } catch( Exception ex )
        {
            throw new RuntimeException( ex );
        } catch( Throwable throwable )
        {
            throw new RuntimeException( throwable );
        }
    }
    
    private void
    readBamFile(ValidationResult result, List<RawReadsFile> files, AtomicBoolean paired )
    {
        readBamOrCramFile(result, files, paired);
    }

    private void
    readBamOrCramFile( ValidationResult result, List<RawReadsFile> files, AtomicBoolean paired )
    {
        BamScanner scanner = new BamScanner()
        {
            @Override protected void
            logProcessedReadNumber( long count )
            {
                ReadsValidator.this.logProcessedReadNumber( count );
            }
        };

        for( RawReadsFile rf : files )
        {
            ValidationResult fileResult;
            if( rf.getReportFile() == null ) {
                fileResult = result.create(new ValidationOrigin("file", rf.getFilename()));
            } else
            {
                fileResult = result.create(rf.getReportFile().toFile(), new ValidationOrigin("file", rf.getFilename()));
            }

            try
            {
                fileResult.add(ValidationMessage.info("Processing file"));
                scanner.readCramFile( fileResult, rf, paired );

            } catch( SAMFormatException | CRAMException ex )
            {
                fileResult.add(ValidationMessage.error(ex));

            } catch( IOException ex )
            {
                throw new RuntimeException( ex );
            }
        }
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
}
