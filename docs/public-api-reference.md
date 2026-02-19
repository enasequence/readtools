# Readtools Public API Reference

This document covers the public interfaces of all classes in the readtools project and describes how they interact.

---

## Class Interactions

### Processing Pipelines

The project implements three main data-flow pipelines and a validation subsystem. All pipelines share common infrastructure for quality normalization, read pairing, and memory management.

### 1. FASTQ-to-BAM Pipeline (`Fastq2Sam`)

```
FASTQ files
  │
  ├─ Utils.detectFastqQualityFormat()        detect Sanger/Illumina/Solexa
  ├─ Utils.getQualityNormalizer()            get matching normalizer
  │
  ├─ ReadReader (one per file)               parse FASTQ records into Read objects
  │     │
  │     └─ ReadWriter chain:
  │          ├─ SingleFastqWriter             (if 1 file) wraps Read → PairedRead
  │          │   OR
  │          ├─ PairedFastqWriter             (if 2 files) pairs reads by key
  │          │   extends AbstractPagedReadWriter
  │          │     └─ spills to disk when memory exceeded
  │          │
  │          └─→ Fastq2BamWriter              downstream writer
  │                ├─ QualityNormalizer.normalize()
  │                ├─ CasavaRead (detect format, extract barcode)
  │                ├─ Utils.replaceUracilBases() (optional)
  │                └─ SAMFileWriter.addAlignment()
  │
  └─ readWriter.cascadeErrors()              flush remaining unpaired reads
```

`Fastq2Sam.convert()` orchestrates this chain. It creates the appropriate `ReadWriter` based on file count, wires it to a `Fastq2BamWriter`, and runs a `MultiFastqConverter` that reads FASTQ records and feeds them through the writer chain. Quality normalization happens inside `Fastq2BamWriter` using a `common.reads.QualityNormalizer` instance selected by `Utils.getQualityNormalizer()`.

The key pairing class is `PairedFastqWriter`, which extends `AbstractPagedReadWriter`. It extracts a read key via `getReadKey()` (regex strips `/1`, `.1`, `:1`, `_1` and Casava tails), buffers reads in a `HashMap<String, List<Read>>`, and assembles `PairedRead` objects when both mates arrive. When memory exceeds the configured limit, `AbstractPagedReadWriter` serializes the map to a GZIP-compressed temp file (`spillMap()`). At `cascadeErrors()` time, spilled maps are reloaded (`fillMap()`) and remaining reads are forwarded as orphans.

### 2. Direct FASTQ Normalization Pipeline (`FastqNormalizer`)

```
FASTQ files
  │
  ├─ Utils.detectFastqQualityFormat()
  ├─ Utils.getQualityNormalizer()
  │
  ├─ FastqReader (htsjdk)                    parse input
  │     │
  │     └─ FastqNormalizer internal logic:
  │          ├─ PairedFastqWriter.getReadKey()     extract read key
  │          ├─ PairedFastqWriter.getPairNumber()  extract pair index
  │          ├─ QualityNormalizer.normalize()       normalize quality bytes
  │          ├─ Utils.replaceUracilBases()          U→T conversion
  │          ├─ CasavaRead.getBaseNameOrNull()      detect Casava format
  │          │
  │          ├─ HashMap<String, List<NormalizedRead>>  in-memory pair buffer
  │          │   └─ spills to disk (ObjectOutputStream + GZIP) on overflow
  │          │
  │          └─ multi-generation spill reassembly
  │                load first spill, stream remaining, merge by key
  │
  └─ AsyncFastqWriter (htsjdk)               write normalized FASTQ output
```

`FastqNormalizer` is a self-contained class that reimplements pairing logic independently of the `ReadWriter` chain. It reuses `PairedFastqWriter`'s static methods (`getReadKey()`, `getPairNumber()`) for read name parsing, and the `common.reads.QualityNormalizer` hierarchy for quality conversion, but does its own buffering and spill-to-disk using Java serialization. This avoids the BAM round-trip that `Fastq2Sam` + `Sam2Fastq` would require.

### 3. BAM/CRAM-to-FASTQ Pipeline (`Sam2Fastq`)

```
SAM/BAM/CRAM file
  │
  ├─ SamReaderFactory
  │   └─ .referenceSource(ENAReferenceSource)   (for CRAM only)
  │
  ├─ SamReader.iterator()                       iterate SAMRecords
  │     │
  │     ├─ skip secondary/supplementary
  │     ├─ skip vendor-QC-failed
  │     ├─ reverse-complement if negative strand
  │     │
  │     └─→ MultiFastqOutputter
  │           ├─ buffers reads, matches pairs by name
  │           ├─ sorts by read name
  │           ├─ writes to 3 output streams: /1, /2, unpaired
  │           └─ overflow: writes excess to temp BAM, re-reads later
  │
  └─ MultiFastqOutputter.finish()               flush remaining
```

`Sam2Fastq` opens the input via htsjdk's `SamReader`. For CRAM files, it provides an `ENAReferenceSource` that resolves reference sequences through a 4-tier cache: in-memory (`SoftReference` map) -> disk (`REF_CACHE` env, falling back to `~/.cache/hts-ref`) -> remote (`REF_PATH` env, default `https://www.ebi.ac.uk/ena/cram/md5/%s`) -> `@SQ:UR` header tag. `ENAReferenceSource` uses `PathPattern` to format MD5-based paths and Spring's `RetryTemplate` for download retries.

`MultiFastqOutputter` collates reads by name and distributes them to paired/unpaired output files. It uses `FastqRead` objects (compact `byte[]`-backed records) that implement `IRead` for `Comparable` sorting.

### 4. Validation Subsystem

```
ReadsManifest (from webin-cli)
  │
  └─ validator.ReadsValidator.validate()
       ├─ convert SubmissionFiles → RawReadsFile list
       │
       └─ v2.validator.ValidatorWrapper(files, format, readLimit)
            │
            ├─ For each file:
            │   ├─ ReadsProviderFactory(file, format)
            │   │   ├─ FASTQ → FastqReadsProvider
            │   │   │             └─ Utils.detectFastqQualityFormat()
            │   │   │             └─ QualityNormalizer for score conversion
            │   │   └─ BAM/CRAM → SamReadsProvider
            │   │                   └─ ENAReferenceSource (for CRAM)
            │   │
            │   ├─ InsdcReadsValidator.validate(factory)
            │   │   └─ checks IUPAC codes, read names, quality thresholds
            │   │
            │   └─ FastqReadsValidator.validate(factory)    (FASTQ only)
            │       ├─ detects read naming style (Casava 1.8 vs generic)
            │       ├─ BloomWrapper for duplicate detection
            │       └─ validates quality score character ranges
            │
            └─ For paired FASTQ (multiple files):
                ├─ PairedFastqReadsValidator per file
                │   └─ BloomWrapper.getCopy() cross-file pairing detection
                │
                └─ calculate pairing percentage (threshold: 20%)
```

The validation chain bridges two package layers. `validator.ReadsValidator` implements the webin-cli `Validator` interface and translates `ReadsManifest` objects into `RawReadsFile` lists. It delegates to `v2.validator.ValidatorWrapper`, which uses the `v2.provider` factory to create format-appropriate `ReadsProvider` instances, then runs `InsdcReadsValidator` (IUPAC codes, quality metrics) and `FastqReadsValidator` (naming, duplicates) against them.

For paired FASTQ validation, `ValidatorWrapper` creates a `BloomWrapper` on the first file, then passes copies (`getCopy()`) to `PairedFastqReadsValidator` instances for each secondary file. The bloom filter detects read names that appear in both files, giving a pairing percentage.

### 5. The `FastqIterativeWriter` Entry Point

```
FastqIterativeWriter.iterator()
  │
  └─ MultiFastqConverterIterator(config, files)
       ├─ SingleFastqWriter or PairedFastqWriter    depending on READ_TYPE
       │   └─ .setWriter(this)                      iterator is the downstream
       │
       ├─ FileCompression.open(file)                decompress input
       │
       └─ MultiFastqConverter(streams, writer)      drives the pipeline
            └─ ReadReader per stream                 FASTQ parsing
```

`FastqIterativeWriter` is the lazy-evaluation entry point used by external consumers. It creates a `MultiFastqConverterIterator` that acts as both `Iterator<PairedRead>` and `ReadWriter<PairedRead, Spot>`. The converter processes one record at a time via `runOnce()`, feeds it through the `ReadWriter` chain, and the iterator's own `write()` method queues assembled `PairedRead` objects. `hasNext()` drives the converter until a record appears in the queue or input is exhausted.

### Shared Infrastructure

**Quality Normalization** has two parallel implementations:
- `common.reads.QualityNormalizer` (interface) with `Standard/Illumina/Solexa` implementations -- used by `Fastq2BamWriter`, `FastqNormalizer`, and `FastqReadsProvider` for byte-level Phred conversion.
- `loader.common.QualityNormalizer` (enum) with `SANGER/SOLEXA/ILLUMINA_1_3/1_5` -- used by `MultiFastqConverter` for string-level quality normalization during initial FASTQ parsing.

`Utils.detectFastqQualityFormat()` auto-detects which format is in use and `Utils.getQualityNormalizer()` returns the matching `common.reads.QualityNormalizer`.

**Read Name Parsing** is centralized in two places:
- `PairedFastqWriter.getReadKey()` / `getPairNumber()` -- used by both `PairedFastqWriter` itself and `FastqNormalizer` for extracting the base read name and pair index from separator formats (`/1`, `.1`, `:1`, `_1`).
- `CasavaRead` -- used by `Fastq2BamWriter`, `FastqNormalizer`, `SingleFastqWriter`, and validators for detecting and parsing Casava 1.8 headers.

**Memory Management** via spill-to-disk uses two strategies:
- `AbstractPagedReadWriter.spillMap()` / `fillMap()` -- Java serialization + GZIP, used by `PairedFastqWriter` in the BAM pipeline.
- `FastqNormalizer`'s internal spill -- similar serialization but with multi-generation reassembly for merging pairs that were split across spill files.

---

## Package: `common.reads`

### `QualityNormalizer` (Interface)
Normalizes quality score byte arrays in-place.

```java
void normalize(byte[] qualities)
```

### `StandardQualityNormalizer` implements `QualityNormalizer`
Converts Sanger FASTQ quality characters (Phred+33) to binary Phred scores via `SAMUtils.fastqToPhred()`.

### `IlluminaQualityNormalizer` implements `QualityNormalizer`
Converts Illumina/Solexa 1.3+ quality characters (Phred+64) to binary Phred scores via `SolexaQualityConverter`.

### `SolexaQualityNormalizer` implements `QualityNormalizer`
Converts Solexa quality characters to binary Phred scores via `SolexaQualityConverter`.

### `CasavaRead`
Static utility for Casava 1.8 read name parsing. Format: `@instrument:run:flowcell:lane:tile:x:y readnum:filter:control:barcode`

```java
static String getBaseNameOrNull(String readName)   // instrument part (group 1), or null
static String getReadIndexOrNull(String readName)   // read number (1, 2, ...), or null
static String getBarcodeOrNull(String readName)     // barcode sequence, or null
static boolean isFiltered(String readName)          // true if filter flag is Y
```

**Constants:**
- `Pattern P_CASAVA_18_NAME` -- pattern for processed read names (no @ prefix)
- `Pattern P_CASAVA_18_RAW_LINE` -- pattern for raw FASTQ lines (with @ prefix)

---

## Package: `cram`

### `CramTools`
CLI entry point for CRAM tools. Registers commands (Sam2Fastq) and dispatches via JCommander.

```java
static void main(String[] args)
```

### `cram.common.Utils`
Reference sequence utilities.

```java
static void checkRefMD5(SAMSequenceDictionary d, ReferenceSequenceFile refFile, boolean checkExistingMD5, boolean failIfMD5Mismatch)
static String calculateMD5String(byte[] data)
static String calculateMD5String(byte[] data, int offset, int len)
static byte[] calculateMD5(byte[] data, int offset, int len)
static int readInto(ByteBuffer buf, InputStream inputStream)
static boolean isValidSequence(byte[] bases, int checkOnlyThisManyBases)
static byte upperCase(byte base)
static byte[] upperCase(byte[] bases)
```

### `InMemoryReferenceSequenceFile` implements `ReferenceSequenceFile`
In-memory reference sequence store with SAMSequenceDictionary.

```java
void addSequence(String name, byte[] bases)
ReferenceSequence getSequence(String name)
SAMSequenceDictionary getSequenceDictionary()
ReferenceSequence getSubsequenceAt(String name, long start, long stop)
ReferenceRegion getRegion(String name, long start, long stop)
boolean isIndexed()          // always true
ReferenceSequence nextSequence()
void reset()
void close()
```

### `PathPattern`
Formats path patterns with hierarchical variable substitution. Example: `"%2s/%3s/%s"` with `"0123456789"` yields `"01/234/56789"`.

```java
PathPattern(String format)
String format(String id)
```

### `ReferenceRegion`
Region of reference bases with alignment coordinates.

```java
// Fields
int index; String name; long alignmentStart; byte[] array;

ReferenceRegion(byte[] bases, int sequenceIndex, String sequenceName, long alignmentStart)
static ReferenceRegion wrap(byte[] bases, int sequenceIndex, String sequenceName)
static ReferenceRegion copyRegion(byte[] bases, int sequenceIndex, String sequenceName, long alignmentStart, long alignmentEnd)
int arrayPosition(long alignmentPosition)
byte base(long alignmentPosition)
byte[] copy(long alignmentStart, int alignmentSpan)
byte[] copySafe(long alignmentStart, int alignmentSpan)
static byte[] copySafe(byte[] ref, long refAlStart, int alignmentStart, int alignmentSpan)
String md5(int alignmentStart, int alignmentSpan)
```

### `ENAReferenceSource` implements `CRAMReferenceSource`
Multi-level caching reference discovery: memory (SoftReferences) -> disk (REF_CACHE) -> remote (REF_PATH, default EBI endpoint) -> @SQ:UR URL fallback.

```java
ENAReferenceSource()
ENAReferenceSource(String... cachePatterns)
synchronized byte[] getReferenceBases(SAMSequenceRecord record, boolean tryNameVariants)
byte[] getReferenceBasesByRegion(SAMSequenceRecord record, int zeroBasedStart, int requestedRegionLength)
synchronized ReferenceRegion getRegion(SAMSequenceRecord record, int start_1based, int endInclusive_1based)
void setLoggerWrapper(LoggerWrapper log)
void clearMemCache()
List<String> getRefPathList()
List<String> getRefCacheList()
int getDownloadTriesBeforeFailing()
void setDownloadTriesBeforeFailing(int n)
static List<String> splitRefPath(String paths)
static RetryTemplate getRetryTemplate(Duration minBackoff, Duration maxBackoff, Double multiplier, int attempts)
```

**Inner interface:** `LoggerWrapper` -- logging abstraction with default no-op methods.

### `CramReferenceException` extends `RuntimeException`
Unchecked exception for CRAM reference errors.

---

## Package: `fastq`

### `IRead` (Interface) extends `Comparable<IRead>`
Marker/ordering interface for FASTQ reads used by the sorter/outputter pipeline.

### `FastqRead` implements `IRead`
Compact byte-array-backed FASTQ record. Stores `@name\nbases\n+\nscores\n` as a single `byte[] data`.

```java
FastqRead(int readLength, byte[] name, boolean appendSegmentIndex, int templateIndex, byte[] bases, byte[] scores)
int compareTo(IRead read)
SAMRecord toSAMRecord(SAMFileHeader header)
```

### `MultiFastqOutputter`
Writes sorted FASTQ reads to paired/unpaired output files with read name deduplication.

```java
MultiFastqOutputter(File[] files, boolean appendSegmentIndex, boolean convertUracil)
MultiFastqOutputter(OutputStream[] streams, boolean appendSegmentIndex, boolean convertUracil)
void write(FastqRead read)
void finish()
long getCount()
```

### `FastqNormalizer`
Direct FASTQ-to-FASTQ normalizer. Quality normalization, read pairing, pair separator canonicalization (`/N`, `.N`, `:N`, `_N` -> `/1`, `/2`), uracil conversion, spill-to-disk for large datasets.

```java
// Paired-end
static PairedNormalizationResult normalizePairedEnd(
    String[] inputFiles, String outputPair1, String outputPair2, String outputSingleton,
    long memLimitBytes, File tmpDir)

// Single-end
static long normalizeSingleEnd(String inputFile, String outputFile)
static SingleNormalizationResult normalizeSingleEndWithStats(String inputFile, String outputFile)
```

**Result classes:**

```java
class PairedNormalizationResult {
    long getPairCount()
    long getOrphanCount()
    long getBaseCount()
}

class SingleNormalizationResult {
    long getReadCount()
    long getBaseCount()
}
```

### `ena.Fastq2Sam`
Converts FASTQ files to BAM via an intermediate sorting step. Uses `MultiFastqOutputter` + `Sam2Fastq` under the hood.

```java
Fastq2Sam()
long convert(String[] inputFiles, String outputFile, String sampleAlias, long memLimitBytes, File tmpDir)
```

### `ena.Fastq2BamWriter`
Writes SAM/BAM records from an iterator of `PairedRead` objects. Handles read pairing flags, quality normalization, Casava read name restoration, and uracil conversion.

```java
Fastq2BamWriter(File outputFile, String sampleAlias)
long write(Iterator<PairedRead> iterator, boolean isConvertUracil) throws IOException
```

---

## Package: `loader.common`

### `FileCompression` (Enum)
Values: `BZ2`, `GZIP`, `GZ`, `ZIP`, `BGZIP`, `BGZ`, `NONE`

```java
InputStream open(String f_name, boolean tar)
InputStream open(File file, boolean tar)
static FileCompression getCompressor(File file)
static InputStream open(File file)
static boolean getUseTar(File file)
```

### `Pair<K, V>` implements `Serializable`
Simple key-value pair with public fields `key` and `value`.

### `QualityNormalizer` (Enum)
Values: `SANGER(33,0,93)`, `SOLEXA(64,-5,93)`, `ILLUMINA_1_3(64,0,93)`, `ILLUMINA_1_5(64,3,93)`, `X(33,0,93)`, `X_2(64,0,93)`, `NONE(33,0,93)`

```java
byte[] normalize(String value) throws QualityNormaizationException
String denormalize(byte[] value) throws QualityNormaizationException
```

### `QualityConverter` (Enum)
Values: `PHRED`, `LOGODDS`

```java
byte[] convertNormalized(QualityConverter to, byte[] what)
```

### `InvalidBaseCharacterException` extends `RuntimeException`
Thrown when invalid characters are found in base sequences.

```java
String getBases()
Collection<Character> getInvalidCharacters()
```

---

## Package: `loader.common.writer`

### `Spot` (Interface)
Base interface for all read/spot data objects.

```java
String getName()
long getBaseCount()
long getSizeBytes()
```

### `ReadWriter<T1 extends Spot, T2 extends Spot>` (Interface)
Chain-of-responsibility interface for processing spots in a pipeline.

```java
void write(T1 spot) throws ReadWriterException
void cascadeErrors() throws ReadWriterException
void setWriter(ReadWriter<T2, ? extends Spot> readWriter)
```

### `AbstractReadWriter<T1, T2>` implements `ReadWriter<T1, T2>`
Base class that groups spots by key, assembles complete collections, and forwards to downstream writer.

```java
AbstractReadWriter()                    // default map size 1M
AbstractReadWriter(int map_size)        // protected

AbstractReadWriter<T1, T2> setVerbose(boolean verbose)
void setWriter(ReadWriter<T2, ? extends Spot> readWriter)
void write(T1 spot)
void cascadeErrors()

// Abstract methods subclasses must implement:
abstract String getKey(T1 spot)
abstract T2 assemble(String key, List<T1> list)
abstract boolean isCollected(List<T1> list)
abstract T2 handleErrors(String key, List<T1> list)

// Overridable hooks:
void append(List<T1> list, T1 spot)
List<T1> newListBucket()
```

### `AbstractPagedReadWriter<T1, T2>` extends `AbstractReadWriter<T1, T2>`
Adds disk-based spill-over when memory limits are exceeded. Serializes maps to GZIP-compressed temp files.

```java
AbstractPagedReadWriter()
AbstractPagedReadWriter(File tmp_root, int spill_page_size, long spill_page_size_bytes, long spill_abandon_limit_bytes)
File spillMap(Map<String, List<T1>> map)
Map<String, List<T1>> fillMap(File file)
void cascadeErrors()   // processes spilled files before cascading
```

### `ReadWriterException` extends `RuntimeException`
Categorized exception with `ErrorType` enum.

```java
ReadWriterException(String value, ErrorType errorType)
ErrorType getErrorType()

enum ErrorType {
    UNEXPECTED_PAIR_NUMBER, SPOT_DUPLICATE, SORTING_ERROR,
    BASES_QUALITIES_LENGTH_MISMATCH, SAM_RECORD_ERROR, INVALID_READ_NAME
}
```

### `ReadWriterMemoryLimitException` extends `RuntimeException`
Thrown when spill memory limits are exceeded.

---

## Package: `loader.fastq`

### `Read` implements `Serializable, Spot`
Single unpaired FASTQ read.

```java
Read(String name, String bases, String qualityScores)
Read(String name, String bases, String qualityScores, String defaultReadIndex)
String getName()
long getBaseCount()
String getBases()
String getQualityScores()
String getDefaultReadIndex()
long getSizeBytes()
```

### `PairedRead` implements `Spot`
A pair of reads (forward/reverse) or a single unpaired read.

```java
// Fields
public final String name;
public final Read forward;
public final Read reverse;

PairedRead(String name, Read read)                    // unpaired
PairedRead(String name, Read forward, Read reverse)   // paired
String getName()
long getBaseCount()
boolean isPaired()
Read getUnpaired()
long getSizeBytes()
```

### `SamRecordWrapper` extends `Read`
Wraps an htsjdk `SAMRecord` as a `Read`, extracting name/bases/qualities.

```java
SamRecordWrapper(SAMRecord samRecord)
SAMRecord getSamRecord()
```

### `PairedFastqWriter` extends `AbstractPagedReadWriter<Read, PairedRead>`
Pairs reads by extracting read keys via regex (`SPLIT_REGEXP` handles `/`, `.`, `:`, `_` separators and Casava format).

```java
PairedFastqWriter(File tmp_root, int spill_page_size, long spill_page_size_bytes, long spill_abandon_limit_bytes)
static String getReadKey(String readname)      // extract base name
static String getPairNumber(String readname)   // extract "1" or "2"
String getKey(Read spot)
PairedRead assemble(String key, List<Read> list)
boolean isCollected(List<Read> list)
PairedRead handleErrors(String key, List<Read> list)
```

### `SingleFastqWriter` implements `ReadWriter<Read, PairedRead>`
Wraps single `Read` objects as unpaired `PairedRead` objects.

```java
void write(Read spot)
void cascadeErrors()
void setWriter(ReadWriter<PairedRead, ?> readWriter)
```

### `FastqIterativeWriter` implements `Iterable<PairedRead>`
Main entry point for FASTQ processing. Configures and iterates over paired reads from input files.

```java
Iterator<PairedRead> iterator()
void setFiles(File[] files)
File[] getFiles()
int getSpillPageSize()
void setSpillPageSize(int spill_page_size)
long getSpill_page_size_bytes()
long getSpill_abandon_limit_bytes()
File getTmpFolder()
void setTmpFolder(File tmp_folder)
READ_TYPE getReadType()                        // SINGLE or PAIRED
void setReadType(READ_TYPE read_type)
QualityNormalizer[] getNormalizers()
void setNormalizers(QualityNormalizer[] normalizers)
Long getReadLimit()
void setReadLimit(Long readLimit)
```

### `MultiFastqConverterIterator` implements `Iterator<PairedRead>, ReadWriter<PairedRead, Spot>`
Iterator that processes FASTQ files via writer pipelines and maintains a queue of assembled `PairedRead` objects.

```java
MultiFastqConverterIterator(File tmp_folder, int spill_page_size, long spill_page_size_bytes,
    long spill_abandon_limit_bytes, READ_TYPE read_type, File[] files,
    QualityNormalizer[] normalizers, Long readLimit)
boolean hasNext()
PairedRead next()
void write(PairedRead spot)
void cascadeErrors()
```

---

## Package: `sam`

### `Sam2Fastq`
Converts SAM/BAM/CRAM to FASTQ. Handles reverse-complement restoration, quality reversal, and read sorting.

```java
// Constants
static final boolean INCLUDE_NON_PF_READS = false
static final boolean INCLUDE_NON_PRIMARY_ALIGNMENTS = false

void convert(File inputFile, File outputFile, long memLimitBytes, File tmpDir)
```

Uses `MultiFastqOutputter` for output.

---

## Package: `sampler`

### `AbstractSqueezerStream`
Abstract base for quality/base compression streams.

### `BaseSqueezer` extends `AbstractSqueezerStream`
Compresses DNA base sequences using dictionary encoding.

### `ScoreSqueezer` extends `AbstractSqueezerStream`
Compresses quality scores using dictionary encoding.

### `bits.BitInputStream` / `bits.BitOutputStream`
Bit-level I/O streams for compressed data.

### `bits.ConstantLengthDataDictionary`
Fixed-width dictionary for bit-level encoding.

### `bits.DictionaryInputStream` / `bits.DictionaryOutputStream`
Dictionary-based compression/decompression streams.

### `bits.BaseDict` / `bits.ScoreDict`
Predefined dictionaries for DNA bases and quality scores.

### `generator.FastqGenerator`
Generates synthetic FASTQ data for testing.

### `intervals.QualInterval` / `intervals.QualScoring`
Quality score interval and scoring utilities for the sampler.

---

## Package: `utils`

### `Utils`
General-purpose utilities for FASTQ processing.

```java
static void replaceUracilBasesInFastq(String inputFastq, String outputFastq) throws IOException
static String replaceUracilBases(String bases)
static FastqQualityFormat detectFastqQualityFormat(String fastqFile1, String fastqFile2)
static QualityNormalizer getQualityNormalizer(FastqQualityFormat qualityType)
static InputStream openFastqInputStream(Path path)   // auto-detects gz/bz2
```

Uracil replacement is case-preserving: `U` -> `T`, `u` -> `t`.

### `SAMReverseReadCheck`
CLI tool to find the first reverse-strand read in a SAM/BAM/CRAM file.

```java
SAMReverseReadCheck(File samFile)
Output check()   // returns first reverse read's name/bases/qualities/flags, or null
static void main(String[] args)
```

### `Sam2FastqStatsComparator`
CLI tool to compare base and quality distributions between two sets of FASTQ files (3 files each: paired1, paired2, unpaired).

```java
static void main(String[] args)   // 6 args: file1a file1b file1c file2a file2b file2c
```

---

## Package: `v2`

### `FileFormat` (Enum)
Values: `BAM`, `CRAM`, `FASTA`, `FASTQ`

### `v2.read.IRead` (Interface)
Read abstraction for the v2 API.

```java
String getName()
String getBases()
String getQualityScores()
```

### `v2.read.FastqRead` implements `v2.read.IRead`
FASTQ read with name, bases, quality scores.

```java
FastqRead(String name, String bases, String qualityScores)
```

### `v2.read.SamRead` implements `v2.read.IRead`
SAM record wrapper.

```java
SamRead(SAMRecord samRecord)
boolean isReadLevelQuality()        // true if quality is "*"
boolean hasQualityControlFlag()     // true if 0x200 flag set
SAMRecord getSamRecord()
```

### `v2.provider.ReadsProvider<T extends IRead>` extends `Iterable<T>, AutoCloseable`
Generic reads source interface.

### `v2.provider.ReadsProviderFactory`
Creates appropriate `ReadsProvider` based on file format.

```java
ReadsProviderFactory(File file, FileFormat format)
ReadsProviderFactory(File file, FileFormat format, boolean normaliseFastqQualityScores)
ReadsProvider<? extends IRead> makeReadsProvider()
```

### `v2.provider.FastqReadsProvider` implements `ReadsProvider<FastqRead>`
Provides FASTQ reads with optional quality normalization.

```java
FastqReadsProvider(File fastqFile) throws ReadsValidationException
FastqReadsProvider(File fastqFile, boolean normaliseQualityScores) throws ReadsValidationException
FastqQualityFormat getQualityFormat()
```

### `v2.provider.SamReadsProvider` implements `ReadsProvider<SamRead>`
Provides SAM/BAM/CRAM reads.

```java
SamReadsProvider(File samFile)
static boolean isCram(File file)   // checks magic number
```

### `v2.validator.ReadsValidationException` extends `Exception`
Validation error with read context.

```java
ReadsValidationException(String errorMessage)
ReadsValidationException(String errorMessage, long readIndex)
ReadsValidationException(String errorMessage, long readIndex, String readName)
String getErrorMessage()   // formatted with read index/name
long getReadIndex()
```

### `v2.validator.ReadsValidator` (Abstract)
Base class for validators.

```java
ReadsValidator(long readCountLimit)
protected abstract boolean validate(ReadsProviderFactory readsProviderFactory)
```

### `v2.validator.InsdcReadsValidator` extends `ReadsValidator`
Validates INSDC standards: IUPAC codes, read name length (max 256), quality thresholds (50% reads must have avg quality >= 30).

```java
InsdcReadsValidator(long readCountLimit)
boolean validate(ReadsProviderFactory readsProviderFactory)
long getReadCount()
long getHighQualityReadCount()
```

### `v2.validator.FastqReadsValidator` extends `ReadsValidator`
Validates FASTQ naming conventions (Casava 1.8 and generic) with bloom-filter duplicate detection.

```java
FastqReadsValidator(long readCountLimit)
boolean validate(ReadsProviderFactory readsProviderFactory)
```

**Enum:** `ReadStyle { FASTQ, CASAVA18 }`

### `v2.validator.PairedFastqReadsValidator` extends `FastqReadsValidator`
Validates paired-end FASTQ files with pairing percentage tracking.

```java
PairedFastqReadsValidator(long readCountLimit, String providerName, BloomWrapper pairingBloomWrapper, Set<String> labels)
long getAddCount()
```

### `v2.validator.ValidatorWrapper`
Orchestrates full validation of one or more files.

```java
ValidatorWrapper(List<File> files, FileFormat format, long readCountLimit)
void run()
boolean isPaired()
List<FileQualityStats> getFileQualityStats()
void validateFastq(File file)
void validateSam(File file)
```

**Inner class:**
```java
class FileQualityStats {
    File getFile()
    long getReadCount()
    long getHighQualityReadCount()
}
```

---

## Package: `validator`

### `ReadsValidator` implements `Validator<ReadsManifest, ReadsValidationResponse>`
Top-level validator for webin-cli integration.

```java
ReadsValidationResponse validate(ReadsManifest manifest)
ReadsValidationResponse validate(File reportFile, ReadsManifest.QualityScore qualityScore,
    SubmissionFiles<FileType> submissionFiles, boolean isQuick)
ReadsValidationResponse validate(ValidationResult result, List<RawReadsFile> files, boolean isQuick)
List<ValidatorWrapper.FileQualityStats> getFileQualityStats()
static RawReadsFile submissionFilesToRawReadsFile(SubmissionFile<FileType> file)
```

Read limits: quick = 100,000; extended = 100,000,000.

---

## Package: `webin.cli.rawreads`

### `RawReadsFile`
File metadata container with type, checksum, quality encoding, and compression info.

**Enums:**
- `Filetype { fastq, bam, cram }`
- `ChecksumMethod { MD5, SHA256 }`
- `QualityScoringSystem { phred, log_odds }`
- `AsciiOffset { FROM33, FROM64 }`
- `QualityEncoding { ascii }`
- `Compression { NONE, GZ, GZIP, BZ2, ZIP }`

Standard getter/setter pairs for `filename`, `filetype`, `checksum`, `checksumMethod`, `qualityScoringSystem`, `asciiOffset`, `qualityEncoding`, `inputDir`, `reportFile`, `compression`.

### `RawReadsException` extends `RuntimeException`
```java
RawReadsException(Throwable ex, String message, String... messages)
RawReadsException(String message, String... messages)
```

### `BloomWrapper`
Bloom filter for duplicate/pairing detection.

```java
BloomWrapper(long expectedReads)
BloomWrapper(long expectedReads, int possibleDuplicatesRetainLimit)
void add(String readName)
long getAddCount()
boolean hasPossibleDuplicates()
Long getPossibleDuplicateCount()
Set<String> getPossibleDuplicates()
boolean contains(String readName)
BloomWrapper getCopy()
Map<String, Set<String>> findAllduplications(String[] read_names, int limit)
```

### `DelegateIterator<T1, T2>` implements `Iterator<T2>`
Abstract type-converting iterator wrapper.

```java
DelegateIterator(Iterator<T1> iterator)
boolean hasNext()
T2 next()
abstract T2 convert(T1 obj)
```

### `InsdcStandardCheckingScanner` implements `ReadWriter<Read, Spot>`
Abstract base for INSDC validation scanners. Validates bases/quality length matching.

```java
InsdcStandardCheckingScanner(int printFreq)
void write(Read read)
protected abstract void logProcessedReadNumber(long cnt)
```

### `FastqReadScanner` extends `InsdcStandardCheckingScanner`
FASTQ-specific scanner with pairing and duplicate detection.

```java
FastqReadScanner(String streamName, Set<String> labels, BloomWrapper pairingBloomWrapper,
    BloomWrapper duplicationsBloomWrapper, int maxLabelSetSize, int printFreq)
void write(Read read)
```

### `FastqScanner` (Abstract)
High-level multi-file FASTQ scanner with pairing analysis.

```java
// Constants
int MAX_LABEL_SET_SIZE = 10
int PAIRING_THRESHOLD = 20

FastqScanner(int expected_size)
FastqScanner(Long readLimit)
FastqScanner(Long readLimit, int expected_size)
boolean getPaired()
void checkFiles(ValidationResult validationResult, RawReadsFile... rawReadsFiles)
protected abstract void logFlushMsg(String message)
protected abstract void logProcessedReadNumber(Long count)
```

### `CramReferenceInfo`
Fetches and validates CRAM reference metadata from the ENA REST service.

```java
ReferenceInfo fetchReferenceMetadata(String md5, PrintStream... pss)
Map<String, Boolean> confirmFileReferences(File file)
```

**Inner class:**
```java
class ReferenceInfo {
    String getId()
    String getMd5()
}
```
