# FASTQ Normalization: Architecture and Data Flow

This document describes the two FASTQ normalization paths in readtools:
the original **pipeline path** (Fastq2Sam -> Sam2Fastq) and the newer
**direct path** (FastqNormalizer). Both produce normalized FASTQ output
with standardized quality scores, read names, and optional uracil
conversion.

## Overview

```
Pipeline path (original):

  input_1.fastq ─┐                                   ┌─> output_1.fastq
                 ├─> Fastq2Sam ─> BAM ─> Sam2Fastq ──┤
  input_2.fastq ─┘                                   ├─> output_2.fastq
                                                     └─> output.fastq (unpaired)

Direct path (new):

  input_1.fastq ─┐                        ┌─> output_1.fastq
                 ├─> FastqNormalizer ─────┤
  input_2.fastq ─┘                        └─> output_2.fastq (pairs + orphans)
```

Both paths perform the same core operations:

1. Quality score detection and normalization (Solexa/Illumina -> Phred+33)
2. Read pairing by key extraction
3. Read name rewriting with optional prefix
4. Pair separator normalization (`.N`, `_N`, `:N` -> `/1`, `/2`)
5. Optional uracil base conversion (U -> T, u -> t)
6. Spill-to-disk for large datasets that exceed memory


## Quality Normalization

Both paths detect the FASTQ quality format using `Utils.detectFastqQualityFormat()`,
which samples reads from the input file(s) and determines the encoding:

| Format    | Offset | ASCII range |
|-----------|--------|-------------|
| Standard (Sanger/Illumina 1.8+) | +33 | 33-126 |
| Illumina (1.3-1.7) | +64 | 64-126 |
| Solexa    | +64 (different formula) | 59-126 |

A format-specific `QualityNormalizer` is selected via
`Utils.getQualityNormalizer()`:

- **StandardQualityNormalizer** -- pass-through (already Phred+33)
- **IlluminaQualityNormalizer** -- subtracts 31 from each byte
- **SolexaQualityNormalizer** -- applies the Solexa-to-Phred conversion formula

After normalization, quality bytes are Phred+33 encoded. The pipeline
path gets this for free via the BAM round-trip (htsjdk handles it).
The direct path calls `normalizer.normalize(qualityBytes)` followed
by `SAMUtils.phredToFastq()`.


## Read Name Handling

### Key Extraction

Read pairing requires extracting a "key" (base name) from each read name.
Both paths use `PairedFastqWriter.getReadKey()`, which tries two strategies
in order:

1. **Casava 1.8 detection** via `CasavaRead.getBaseNameOrNull()`:
   Matches the pattern `<instrument> <index>:<filter>:<control>:<barcode>`.
   Example: `A00500:310:HG3JFDRXY:1:2101:1000:2000 1:N:0:ATCACG`
   -> key = `A00500:310:HG3JFDRXY:1:2101:1000:2000`

2. **Generic SPLIT_REGEXP** `^(.*)(?:[.|:|/|_])([0-9]+)$`:
   Splits on the last `.`, `:`, `/`, or `_` followed by digits.
   A `CASAVA_LIKE_EXCLUDE_REGEXP` guard prevents Casava-like names
   (7 colon-separated fields like `A00953:544:HMTFHDSX3:2:1101:6768:1`)
   from being incorrectly split by this regex.

   Examples:
   - `MYREAD/1`  -> key = `MYREAD`, index = `1`
   - `MYREAD.2`  -> key = `MYREAD`, index = `2`
   - `MYREAD_1`  -> key = `MYREAD`, index = `1`
   - `MYREAD:1`  -> key = `MYREAD`, index = `1`
   - `MYREAD.10` -> key = `MYREAD`, index = `10`

### Output Read Name Format

Both paths normalize the output read name to:

```
{prefix}.{counter} {baseName}/{pairIndex}
```

Or without prefix:

```
{baseName}/{pairIndex}
```

Where:
- `prefix` is typically a run accession (e.g. `ERR12345`)
- `counter` is a sequential number starting at 1
- `baseName` is the key extracted above
- `pairIndex` is `1` or `2` (always, regardless of original separator or index value)

**All separator styles are normalized to `/`.**
Multi-digit pair indices (e.g. `.10`, `.20`) are normalized to `/1`, `/2`
based on the order of first appearance.

### Casava Read Names

Casava 1.8 names receive special treatment:

- **Pipeline path**: The BAM round-trip strips the Casava metadata tail
  (the ` 1:N:0:ATCACG` part after the space). Output becomes
  `ERR.1 A00500:310:HG3JFDRXY:1:2101:1000:2000/1`.

- **Direct path**: FastqNormalizer preserves the full Casava name as-is
  (no `/1`/`/2` appended, since the pair number is already in the Casava
  tail). Output becomes
  `ERR.1 A00500:310:HG3JFDRXY:1:2101:1000:2000 1:N:0:ATCACG`.

This is an intentional behavioral difference.


## Pipeline Path: Fastq2Sam -> Sam2Fastq

### Fastq2Sam

Entry point: `Fastq2Sam.create(Params)`.

```
Fastq2Sam.Params:
  files              -- 1 or 2 input FASTQ paths
  compression        -- input compression (NONE, GZ, BZ2, ZIP, BGZ)
  data_file          -- output BAM path
  sample_name        -- SAM @RG SM tag (required)
  tmp_root           -- temp directory for spill files
  convertUracil      -- U/u -> T/t conversion
  spill_page_size    -- max entries in memory (default 4,500,000)
  spill_page_size_bytes -- max memory bytes (default 4 GB)
  spill_abandon_limit_bytes -- max spill disk usage (default 10 GB)
```

Data flow for paired-end:

```
InputStreams -> MultiFastqConverter -> PairedFastqWriter -> Fastq2BamWriter -> BAM
```

- **MultiFastqConverter**: Reads FASTQ records from one or two input
  streams in round-robin. Passes each `Read` to the writer chain.
  Tracks read count and base count. Stops when the shorter file hits
  EOF (via `ConverterEOFException`).

- **PairedFastqWriter** (extends `AbstractPagedReadWriter`): Buffers
  reads in a `Map<String, List<Read>>` keyed by `getReadKey()`. Each
  entry holds a 2-element list (slot 0 for first mate, slot 1 for
  second). When both slots are filled, the pair is assembled into a
  `PairedRead` and passed to the downstream writer.

- **Fastq2BamWriter**: Converts `PairedRead` / `Read` objects into
  `SAMRecord`s with proper flags (paired, first/second of pair, mate
  unmapped for orphans) and writes them to a BAM file via htsjdk.

### Sam2Fastq

Entry point: `Sam2Fastq.create(Params)`.

```
Sam2Fastq.Params:
  samFile        -- input BAM/CRAM path
  fastqBaseName  -- output path prefix
  prefix         -- read name prefix (e.g. run accession)
  nofStreams     -- number of output streams (typically 3)
  gzip           -- gzip output
  reverse        -- reverse-complement negative strand reads
```

Output files (with `nofStreams=3`):
- `{baseName}.fastq` -- stream 0: unpaired / orphan reads
- `{baseName}_1.fastq` -- stream 1: first of pair
- `{baseName}_2.fastq` -- stream 2: second of pair

Data flow:

```
BAM -> SamReader -> CollatingDumper -> MultiFastqOutputter -> FASTQ files
```

- **MultiFastqOutputter**: Re-pairs reads using an in-memory cache
  (`TreeMap`). When a mate collision is found, both reads are written
  to their respective output streams. Unpaired reads go to stream 0.
  Overflow records (cache evictions) are written to a temporary BAM
  and re-processed.

  The `getSegmentIndexInTemplate(flags)` method determines the output
  stream: `0` for unpaired, `1` for first-of-pair (flag 0x40), `2`
  for second-of-pair.

### Side Effects of BAM Round-Trip

The BAM intermediate format causes some transformations:

- **Base uppercasing**: BAM stores bases as 4-bit encoded values; on
  read-back they come out uppercase. Lowercase bases in the input are
  lost (including soft-masking information).
- **Casava tail stripping**: BAM read names don't preserve the space-
  separated Casava metadata tail. It is lost in the round-trip.
- **Quality re-encoding**: Quality scores pass through htsjdk's Phred
  encoding, ensuring Phred+33 output.


## Direct Path: FastqNormalizer

Entry point (paired-end): `FastqNormalizer.normalizePairedEnd(...)`.
Entry point (single-end): `FastqNormalizer.normalizeSingleEnd(...)`.

```
FastqNormalizer.normalizePairedEnd params:
  inputFastq1, inputFastq2   -- input FASTQ paths
  outputFastq1, outputFastq2 -- output FASTQ paths
  prefix                     -- read name prefix (nullable)
  convertUracil              -- U/u -> T/t conversion
  tempDir                    -- temp directory for spill files
  spillPageSize              -- max entries in memory (default 100,000)
  spillPageSizeBytes         -- max memory bytes (default 4 GB)
  spillAbandonLimitBytes     -- max spill disk usage (default 10 GB)
```

Returns `PairedNormalizationResult` with `pairCount`, `orphanCount`,
`baseCount`, `totalReadCount`.

### Data Flow (Paired-End)

```
input_1.fastq ─┐
               ├─> processInputFiles() -> pairMap (in memory)
input_2.fastq ─┘         |
                          ├─ if fits in memory -> writeFromMemory()
                          └─ if spilled -> processSpillFiles() -> writeFromMemory()
                                                    |
                                           output_1.fastq, output_2.fastq
```

The `PairedNormalizer` inner class manages the process:

1. **processInputFiles()**: Reads both FASTQ files round-robin using
   iterators. Each record is processed:
   - Bases: optional uracil conversion via `Utils.replaceUracilBases()`
   - Quality: normalize via `QualityNormalizer` + `SAMUtils.phredToFastq()`
   - Key: extract via `PairedFastqWriter.getReadKey()`
   - Pair number: extract via `PairedFastqWriter.getPairNumber()`
   - Store as `NormalizedRead` in `pairMap<String, List<NormalizedRead>>`
     (same 2-slot structure as PairedFastqWriter)
   - If memory limits exceeded, spill to disk via `spillToDisk()`

2. **writeFromMemory()**: Sorts keys lexicographically (matching BAM
   queryname sort order). For each key:
   - Both slots filled -> `writePair()` to output_1 and output_2
   - One slot null -> `writeOrphan()` to output_1

3. **processSpillFiles()** (if spilling occurred): Multi-generation
   merge following the `AbstractPagedReadWriter.cascadeErrors()` pattern:
   - Load one spill file into pairMap as the base
   - Stream remaining files entry-by-entry against the base
   - Matches merge into pairMap; non-matches re-spill to next generation
   - After each generation, write completed pairs/orphans
   - Repeat until no more generation files remain

### Spill File Format

Spill files are gzipped Java `ObjectOutputStream` streams containing
`Pair<String, List<NormalizedRead>>` entries. `NormalizedRead` implements
`Serializable`. Files are stored in `tempDir` and deleted after processing.

### Differences from Pipeline Path

| Aspect | Pipeline (Fastq2Sam -> Sam2Fastq) | Direct (FastqNormalizer) |
|--------|-----------------------------------|--------------------------|
| Intermediate format | BAM file on disk | In-memory map + spill files |
| Base casing | Uppercased by BAM | Preserved as-is |
| Casava names | Tail stripped by BAM | Preserved in full |
| EOF handling | Stops when shorter file ends | Reads both files to completion |
| Default spill page size | 4,500,000 entries | 100,000 entries |
| Output files | 3 files (paired_1, paired_2, unpaired) | 2 files (orphans go to file 1) |
| Performance | BAM encode + decode overhead | Direct FASTQ-to-FASTQ, no BAM |


## Uracil Conversion

Both paths use `Utils.replaceUracilBases()` for U -> T conversion.
The conversion preserves case: uppercase `U` becomes `T`, lowercase
`u` becomes `t`. This preserves any soft-masking information encoded
in base casing.

The conversion is applied per-read during ingestion (before storage
in the pair map or BAM).


## Spill-to-Disk Mechanism

Both paths use a similar spill mechanism for handling datasets too large
to fit in memory. The mechanism is controlled by three parameters:

- **spill_page_size**: Maximum number of map entries before spilling.
- **spill_page_size_bytes**: Maximum estimated memory usage before spilling.
- **spill_abandon_limit_bytes**: Maximum total bytes written to spill files.
  Exceeding this throws `ReadWriterMemoryLimitException`.

### Multi-Generation Reassembly

When mates end up in different spill files (e.g., read A/1 in file 1
and read A/2 in file 3), multi-generation reassembly handles this:

```
Generation 0: [spill_0, spill_1, spill_2]
  Load spill_0 into memory
  Stream spill_1: matches merge, non-matches -> spill_3 (generation 1)
  Stream spill_2: matches merge, non-matches -> spill_4 (generation 1)
  Write pairs/orphans from memory

Generation 1: [spill_3, spill_4]
  Load spill_3 into memory
  Stream spill_4: matches merge, non-matches -> spill_5 (generation 2)
  Write pairs/orphans from memory

...repeat until no new generation files...
```

Within each generation, output is sorted lexicographically by read key.
Across generations, the output is sorted per-generation but not globally.
This matches the existing pipeline behavior.


## API Quick Reference

### FastqNormalizer

```java
// Single-end (returns read count only, backward-compatible)
long count = FastqNormalizer.normalizeSingleEnd(
    inputFastq, outputFastq, prefix, convertUracil);

// Single-end with full stats
SingleNormalizationResult result = FastqNormalizer.normalizeSingleEndWithStats(
    inputFastq, outputFastq, prefix, convertUracil);
result.getReadCount();
result.getBaseCount();

// Paired-end (default thresholds)
PairedNormalizationResult result = FastqNormalizer.normalizePairedEnd(
    inputFastq1, inputFastq2, outputFastq1, outputFastq2,
    prefix, convertUracil, tempDir);

// Paired-end (custom thresholds)
PairedNormalizationResult result = FastqNormalizer.normalizePairedEnd(
    inputFastq1, inputFastq2, outputFastq1, outputFastq2,
    prefix, convertUracil, tempDir,
    spillPageSize, spillPageSizeBytes, spillAbandonLimitBytes);

result.getPairCount();
result.getOrphanCount();
result.getTotalReadCount();  // pairCount * 2 + orphanCount
result.getBaseCount();
```

### Fastq2Sam

```java
Fastq2Sam.Params params = new Fastq2Sam.Params();
params.files = Arrays.asList("input_1.fastq", "input_2.fastq");
params.data_file = "output.bam";
params.sample_name = "SM-001";
params.tmp_root = "/tmp";
params.convertUracil = false;

Fastq2Sam f2s = new Fastq2Sam();
f2s.create(params);

f2s.getTotalReadCount();
f2s.getTotalBaseCount();
```

### Sam2Fastq

```java
Sam2Fastq.Params params = new Sam2Fastq.Params();
params.samFile = new File("output.bam");
params.fastqBaseName = "/output/result";  // produces result_1.fastq, result_2.fastq, result.fastq
params.prefix = "ERR12345";
params.nofStreams = 3;

Sam2Fastq s2f = new Sam2Fastq();
s2f.create(params);
```
