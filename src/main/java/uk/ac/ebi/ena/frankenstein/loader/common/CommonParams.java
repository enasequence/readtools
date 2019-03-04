package uk.ac.ebi.ena.frankenstein.loader.common;

import java.util.List;

import com.beust.jcommander.Parameter;

public class 
CommonParams
{
    final public static int PARAM_ERROR = 2;
    final public static int OK_CODE     = 0;
    final public static int FAILURE     = 1;
    
    
    @Parameter( names = { "-h", "--help" } )
    public boolean help = false;
    
    @Parameter( names = { "-f", "--file" }, description = "files to be loaded, repeat option and parameter in case of more than one files (NB: order matters!)" )
    public List<String> files;

    @Parameter( names = { "-c", "--compression"}, required = false, description = "compression (applied for all input files), supported values: BZ2, GZIP or GZ, ZIP, BGZIP or BGZ and NONE" )
    public String compression = FileCompression.NONE.name(); 
    
    @Parameter( names = { "-o", "--output-data-folder" }, description = "Output folder for data base" )
    public String data_folder_name = "data.tmp";
    
    @Parameter( names = { "-q", "--quality-type" }, description = "types: SANGER, SOLEXA, ILLUMINA_1_3, ILLUMINA_1_5, if set - overrides values of -qs and -qe" )
    public String quality_type = null;
    
    @Parameter( names = { "-qe", "--quality-encoding" }, description = "X|X_2" )
    public String quality_encoding = QualityNormalizer.X.toString();

    @Parameter( names = { "-qs", "--quality-scoring" }, description = "PHRED|LOGODDS" )
    public String quality_scoring = QualityConverter.PHRED.toString();

    @Parameter( names = { "-sl", "--spot-length" }, description = "expected spot length, -1 means variable length" )
    public int spot_length = -1;
    
    @Parameter( names = { "-v", "--verbose" }, description = "Verbose" )
    public boolean verbose = false;
    
    @Parameter( names = { "-md5" }, arity=1, description = "Use md5 in output database. Use -md5 false to switch it off" )
    public boolean use_md5 = true;
    
    @Parameter( names = { "-tar" }, description = "Set this flag if input files are tarred" )
    public boolean use_tar = false;

    @Parameter( names = { "-sps", "-spill-page-size" }, description = "Spill page size, depends on maximum of avaliable memory and size of unassembled record pool" )
    public int spill_page_size = 7000000;

    @Parameter( names = { "-tmp", "--tmp-root" }, description = "Folder to store temporary output in case of paired reads assembly. Should be large enough" )
    public String tmp_root = ".";
    
    
    public String
    toString()
    {
        return String.format( "CommonParams:\nfiles: %s\ncompression: %s\ndata_folder_name: %s\nquality_type: %s\nquality_encoding: %s\nquality_scoring: %s", 
                              files, 
                              compression,
                              data_folder_name,
                              quality_type,
                              quality_encoding,
                              quality_scoring ); 
    }
    
}

