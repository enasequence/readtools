# FastqNormalizer - Direct FASTQ Quality Score Normalization

## Overview

`FastqNormalizer` provides direct FASTQ-to-FASTQ conversion with quality score normalization, avoiding the overhead of BAM intermediate format used in the Fastq2Sam → Sam2Fastq pipeline.

This utility is designed for workflows that need normalized FASTQ output without requiring BAM files, significantly reducing processing time and disk I/O for large datasets.

## Features

- **Quality Score Normalization**: Auto-detects input format (Illumina 1.3+, Solexa, Sanger) and normalizes to Sanger/Phred+33
- **Uracil Base Conversion**: Optional U→T base conversion for RNA sequences
- **Run ID Prefixes**: Optional run ID prefixes matching Sam2Fastq output format
- **Paired-End Support**: Full read pairing, sorting, and validation for paired FASTQ files
- **Memory Management**: Configurable spill-to-disk with abandonment thresholds for large datasets
- **Format Support**: Handles gzip, bzip2, and plain FASTQ files transparently
- **Casava 1.8 Support**: Properly handles Casava 1.8 read name format

## API

### Single-End Normalization

```java
public static long normalizeSingleEnd(
    String inputFastq,
    String outputFastq,
    String prefix,
    boolean convertUracil) throws IOException
```

**Parameters:**
- `inputFastq` - Path to input FASTQ file (gz/bz2/plain auto-detected)
- `outputFastq` - Path to output FASTQ file (extension determines compression)
- `prefix` - Optional run ID prefix (nullable). When set, read names become `{prefix}.{counter} {originalName}`
- `convertUracil` - If true, converts U bases to T

**Returns:** Number of reads written

**Example:**
```java
long count = FastqNormalizer.normalizeSingleEnd(
    "input.fastq.gz",
    "output.fastq.gz",
    "SRR123456",
    true  // convert U to T
);
System.out.println("Wrote " + count + " reads");
```

### Paired-End Normalization (Default Thresholds)

```java
public static long normalizePairedEnd(
    String inputFastq1,
    String inputFastq2,
    String outputFastq1,
    String outputFastq2,
    String prefix,
    boolean convertUracil,
    File tempDir) throws IOException
```

**Parameters:**
- `inputFastq1` - Path to first mate FASTQ file
- `inputFastq2` - Path to second mate FASTQ file
- `outputFastq1` - Path to output first mate FASTQ file
- `outputFastq2` - Path to output second mate FASTQ file
- `prefix` - Optional run ID prefix (nullable)
- `convertUracil` - If true, converts U bases to T
- `tempDir` - Directory for temporary spill files

**Returns:** Number of pairs written

**Default Thresholds:**
- **Spill page size**: 100,000 reads in memory
- **Spill page size bytes**: 4 GB memory threshold
- **Spill abandon limit**: 10 GB disk limit

**Example:**
```java
long pairs = FastqNormalizer.normalizePairedEnd(
    "input_1.fastq.gz",
    "input_2.fastq.gz",
    "output_1.fastq.gz",
    "output_2.fastq.gz",
    "ERR654321",
    false,  // no uracil conversion
    new File("/tmp")
);
System.out.println("Wrote " + pairs + " pairs");
```

### Paired-End Normalization (Custom Thresholds)

```java
public static long normalizePairedEnd(
    String inputFastq1,
    String inputFastq2,
    String outputFastq1,
    String outputFastq2,
    String prefix,
    boolean convertUracil,
    File tempDir,
    int spillPageSize,
    long spillPageSizeBytes,
    long spillAbandonLimitBytes) throws IOException
```

**Additional Parameters:**
- `spillPageSize` - Maximum number of reads to keep in memory before spilling
- `spillPageSizeBytes` - Maximum memory usage in bytes before spilling
- `spillAbandonLimitBytes` - Maximum total spilled bytes before aborting

**Example:**
```java
long pairs = FastqNormalizer.normalizePairedEnd(
    "input_1.fastq.gz",
    "input_2.fastq.gz",
    "output_1.fastq.gz",
    "output_2.fastq.gz",
    "RUN001",
    true,
    new File("/tmp"),
    50_000,                        // 50K reads in memory
    2L * 1024L * 1024L * 1024L,   // 2 GB memory threshold
    5L * 1024L * 1024L * 1024L    // 5 GB disk abandon limit
);
```

## Quality Score Formats

The normalizer automatically detects and converts the following quality score formats to Sanger (Phred+33):

| Format | ASCII Range | Phred Range | Description |
|--------|-------------|-------------|-------------|
| Sanger | 33-126 | 0-93 | Standard FASTQ format (Phred+33) |
| Illumina 1.3+ | 64-126 | 0-62 | Illumina 1.3-1.7 (Phred+64) |
| Solexa | 59-126 | -5-62 | Early Solexa/Illumina (Solexa+64) |

All output FASTQ files use **Sanger (Phred+33)** encoding.

## Read Name Format

### Without Prefix
Input: `@INSTRUMENT:RUN:FLOWCELL:LANE:TILE:X:Y 1:N:0:BARCODE`  
Output: `@INSTRUMENT:RUN:FLOWCELL:LANE:TILE:X:Y/1`

### With Prefix (e.g., "SRR123456")
Input: `@INSTRUMENT:RUN:FLOWCELL:LANE:TILE:X:Y 1:N:0:BARCODE`  
Output: `@SRR123456.1 INSTRUMENT:RUN:FLOWCELL:LANE:TILE:X:Y/1`

The format matches Sam2Fastq output: `{prefix}.{counter} {originalName}/1`

## Paired-End Processing

### Read Pairing Logic

The paired-end normalizer implements full read pairing and sorting:

1. **Read Key Extraction**: Extracts base name from read names (strips `/1`, `/2` suffixes)
2. **Pairing**: Matches reads from both files by their base name
3. **Sorting**: Ensures first-of-pair comes before second-of-pair
4. **Validation**: Detects mismatched file lengths and duplicate reads

### Memory Management

When processing large or out-of-order paired FASTQ files:

1. **In-Memory Buffering**: Reads are buffered in memory up to `spillPageSize` or `spillPageSizeBytes`
2. **Spill to Disk**: When thresholds exceeded, data is serialized to temporary gzip files
3. **Multi-Pass Processing**: Spilled data is processed in multiple passes to reassemble pairs
4. **Abandonment**: If total spilled data exceeds `spillAbandonLimitBytes`, processing aborts with `ReadWriterMemoryLimitException`

### When to Use Custom Thresholds

**Default thresholds work well for:**
- Small to medium datasets (< 100M reads)
- Files with reads in matching order
- Systems with 8+ GB RAM

**Use custom thresholds for:**
- Very large datasets (100M+ reads)
- Severely out-of-order paired files
- Memory-constrained environments
- Systems with fast SSD storage (can use higher disk limits)

## Error Handling

### Common Exceptions

**`IOException`**
- File not found
- Permission denied
- Disk full
- Mismatched paired file lengths

**`ReadWriterMemoryLimitException`**
- Total spilled data exceeds abandonment threshold
- Indicates severely out-of-order paired files
- Solution: Increase `spillAbandonLimitBytes` or pre-sort input files

**`ReadWriterException`**
- Invalid read name format
- Duplicate read in same file
- Unexpected pair number

## Performance Considerations

### Single-End
- **Memory**: Minimal, processes reads sequentially
- **Disk I/O**: 1 read + 1 write per record
- **Speed**: ~1-2M reads/second

### Paired-End (In-Order)
- **Memory**: Buffered up to spillPageSize (default 100K reads ≈ 100-200 MB)
- **Disk I/O**: 2 reads + 2 writes per pair
- **Speed**: ~500K-1M pairs/second

### Paired-End (Out-of-Order with Spilling)
- **Memory**: Capped at spillPageSizeBytes (default 4 GB)
- **Disk I/O**: Additional temporary file I/O proportional to disorder
- **Speed**: ~100K-500K pairs/second (depends on spill frequency)

## Comparison with Fastq2Sam → Sam2Fastq

| Feature | FastqNormalizer | Fastq2Sam → Sam2Fastq |
|---------|----------------|----------------------|
| Quality normalization | ✓ | ✓ |
| Uracil conversion | ✓ | ✓ |
| Read name prefix | ✓ | ✓ |
| Paired-end support | ✓ | ✓ |
| BAM output | ✗ | ✓ |
| BAM encoding overhead | None | ~30-50% slower |
| Disk usage | Input + Output | Input + BAM + Output |
| Memory efficiency | Higher | Lower |
| Use case | FASTQ normalization | Full SAM/BAM workflow |

**Use FastqNormalizer when:**
- You only need normalized FASTQ output
- Processing large datasets where speed matters
- Disk space is limited
- BAM intermediate is not required

**Use Fastq2Sam → Sam2Fastq when:**
- You need BAM output for other tools
- BAM metadata preservation is required (flags, tags)
- Full SAM/BAM pipeline integration is needed

## Examples

### Example 1: Simple Quality Normalization
```java
// Normalize a single FASTQ file from Illumina to Sanger
FastqNormalizer.normalizeSingleEnd(
    "illumina_1.3.fastq.gz",
    "sanger.fastq.gz",
    null,   // no prefix
    false   // no uracil conversion
);
```

### Example 2: RNA-Seq with Uracil Conversion
```java
// Convert RNA reads (U bases) to DNA (T bases) and normalize
FastqNormalizer.normalizeSingleEnd(
    "rna_seq.fastq.gz",
    "dna_normalized.fastq.gz",
    null,
    true    // convert U to T
);
```

### Example 3: Paired-End with Run ID
```java
// Normalize paired-end reads with run ID prefix
FastqNormalizer.normalizePairedEnd(
    "sample_R1.fastq.gz",
    "sample_R2.fastq.gz",
    "ERR123456_1.fastq.gz",
    "ERR123456_2.fastq.gz",
    "ERR123456",  // run ID prefix
    false,
    new File("/tmp")
);
```

### Example 4: Large Dataset with Custom Thresholds
```java
// Process 500M paired reads with higher memory limits
FastqNormalizer.normalizePairedEnd(
    "large_R1.fastq.gz",
    "large_R2.fastq.gz",
    "normalized_1.fastq.gz",
    "normalized_2.fastq.gz",
    "SRR999999",
    false,
    new File("/scratch"),
    500_000,                          // 500K reads in memory
    8L * 1024L * 1024L * 1024L,      // 8 GB memory threshold
    50L * 1024L * 1024L * 1024L      // 50 GB disk limit
);
```

## Integration with webin-read-stages

The FastqNormalizer can be integrated into the webin-read-stages pipeline as an alternative to Fastq2Sam when BAM output is not required:

```java
// In Loading.java
if (skipBAMOutput) {
    // Use FastqNormalizer for direct FASTQ→FASTQ conversion
    FastqNormalizer.normalizePairedEnd(
        inputFastq1,
        inputFastq2,
        outputFastq1,
        outputFastq2,
        runId,
        true,  // convertUracil
        eraConfig.getTempDir().toFile()
    );
} else {
    // Use Fastq2Sam for FASTQ→BAM→FASTQ pipeline
    Fastq2Sam.create(params);
}
```

## Testing

Comprehensive test coverage includes:
- Single-end: Sanger passthrough, quality normalization, prefix handling, uracil conversion
- Paired-end: In-order reads, out-of-order reads, prefix handling, Casava format, mismatched lengths
- Real file testing with actual test resources

Run tests:
```bash
./gradlew test --tests FastqNormalizerTest
```

## See Also

- `uk.ac.ebi.ena.readtools.fastq.ena.Fastq2Sam` - FASTQ to BAM conversion
- `uk.ac.ebi.ena.readtools.sam.Sam2Fastq` - BAM to FASTQ conversion
- `uk.ac.ebi.ena.readtools.utils.Utils` - Quality format detection and normalization utilities
- `uk.ac.ebi.ena.readtools.loader.fastq.PairedFastqWriter` - Paired FASTQ read pairing logic
