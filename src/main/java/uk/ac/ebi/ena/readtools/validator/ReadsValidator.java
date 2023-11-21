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
package uk.ac.ebi.ena.readtools.validator;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.time.Duration;
import java.util.List;
import java.util.Map;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.stream.Collectors;

import htsjdk.samtools.SAMFormatException;
import htsjdk.samtools.cram.CRAMException;

import uk.ac.ebi.ena.readtools.webin.cli.rawreads.SamScanner;
import uk.ac.ebi.ena.readtools.webin.cli.rawreads.FastqScanner;
import uk.ac.ebi.ena.readtools.webin.cli.rawreads.RawReadsFile;
import uk.ac.ebi.ena.readtools.webin.cli.rawreads.RawReadsFile.Filetype;
import uk.ac.ebi.ena.readtools.webin.cli.rawreads.refs.CramReferenceInfo;
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
    private static final Duration DEFAULT_QUICK_RUN_DURATION = Duration.ofMinutes(5);

    private static final long QUICK_READ_LIMIT = 100_000;
    private static final long EXTENDED_READ_LIMIT = 100_000_000;

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
        if (files != null && files.size() > 0) {
            Filetype fileType = files.get(0).getFiletype();

            if( Filetype.fastq.equals( fileType ) )
            {
                readFastqFile( result, files, paired );
            } else if( Filetype.bam.equals( fileType ) )
            {
                readBamFile( result, files, paired );
            } else if( Filetype.cram.equals( fileType ) )
            {
                readCramFile( result, files, paired );
            } else
            {
                throw new RuntimeException( "Unsupported file type: " + fileType.name() );
            }
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
            FastqScanner fs = new FastqScanner(manifest.isQuick() ? QUICK_READ_LIMIT : EXTENDED_READ_LIMIT) {
                @Override
                protected void logFlushMsg( String msg ) {}

                @Override
                protected void logProcessedReadNumber(Long count) {
                    ;
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
        SamScanner scanner = new SamScanner(
                manifest.isQuick() ? QUICK_READ_LIMIT : EXTENDED_READ_LIMIT,
                manifest.isQuick() ? DEFAULT_QUICK_RUN_DURATION : null)
        {
            @Override
            protected void logProcessedReadNumber( long count ) {}
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
                if (rf.getFiletype().equals(Filetype.cram)) {
                    // Cram references have already been confirmed at this point so set the relevant flag to false.
                    scanner.readCramFile( fileResult, rf, paired, false );
                } else {
                    scanner.readBamFile( fileResult, rf, paired );
                }

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
}
