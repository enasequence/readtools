# readtools `validate` CLI

A small standalone command-line wrapper around readtools' read validation, so you can run
**the same client-side FASTQ / BAM / CRAM checks that ENA's `webin-cli` performs** on files you
receive — without the manifest, credentials, or upload machinery.

It is a thin entry point over `uk.ac.ebi.ena.readtools.v2.validator.ValidatorWrapper` (the exact
class `webin-cli` drives via readtools). Source: `src/main/java/uk/ac/ebi/ena/readtools/v2/cli/ValidateCli.java`.

> **Scope:** this reproduces **client-side** validation only. It does not perform ENA's
> server-side / post-submission checks (e.g. technical-read screening, or the average-quality
> rule that is disabled client-side). A file that passes here can still be flagged by ENA after
> upload. Treat it as a fast local pre-flight gate, not a guarantee of acceptance.

---

## Requirements

- **To run the built jar: any JDK ≥ 17.** The jar is compiled to Java 17 bytecode
  (`targetCompatibility = VERSION_17`), which is forward-compatible, so it runs unchanged on
  JDK 17, 21, 25, … (verified on 17 and 21). No newer-JDK-specific build is needed.
- **To build from source: JDK 17–21.** readtools bundles Gradle 7.2; building has been verified
  on JDK 17 and 21. For newer JDKs, bump the Gradle wrapper (`gradle/wrapper/gradle-wrapper.properties`)
  to a version that supports them — the emitted bytecode stays Java 17 regardless of the building
  JDK, so the result still runs on any JDK ≥ 17.
- Internet access to Maven Central for the first build (36 of 37 dependencies resolve from
  Central and are cached afterwards).

## Dependencies — what's not on Maven Central

readtools' entire dependency tree is on Maven Central **except one artifact**:

```
uk.ac.ebi.ena.webin-cli:webin-cli-validator
```

This is a small Apache-2.0 leaf module (interfaces/types only; its own dependencies are all on
Central). It is **not published to Maven Central** and the EBI registry that hosts it is
auth-gated, so you build it from source once and install it into your local Maven repo. Because
the source is Apache-2.0 you can freely clone/mirror it as a backup.

---

## Build

### 1. Build the one off-Central dependency (`webin-cli-validator`) into `~/.m2`

```bash
git clone https://github.com/enasequence/webin-cli-validator.git
cd webin-cli-validator
# APP_VERSION supplies the module version; the -P flag satisfies a property the
# publish block references. Pin an explicit 2.x version for reproducibility.
APP_VERSION=2.15.1 ./gradlew publishToMavenLocal -Pgitlab_private_token=x
cd ..
```

This installs `uk.ac.ebi.ena.webin-cli:webin-cli-validator:2.15.1` into `~/.m2/repository`, where
readtools' build (which lists `mavenLocal()` first) will find it.

> **Tip — pin the version.** readtools declares `webin-cli-validator:2.+` (a floating range),
> which can resolve unpredictably. In your fork, pin it to the exact version you built above
> (e.g. `2.15.1`) so builds are reproducible and don't depend on any remote registry.
>
> **Alternative — no `~/.m2` step:** use a Gradle composite build. Add
> `includeBuild('../webin-cli-validator')` to a `settings.gradle` in this repo and Gradle will
> compile the validator from source on demand — no publish, no version pin, no registry.

### 2. Build the CLI fat jar

```bash
git clone https://github.com/<your-org>/readtools.git   # your fork
cd readtools
./gradlew shadowJar
```

Output (the `shadow` plugin is already configured in `build.gradle`):

```
build/libs/readtools-<version>-all.jar     # ~19 MB, self-contained
```

The jar's basename follows the project directory name; for a fork checked out as `readtools`
it is `readtools-2.15.1-all.jar`.

---

## Usage

The fat jar is directly executable (its manifest sets `Main-Class`):

```bash
java -jar build/libs/readtools-2.15.1-all.jar --format FASTQ [--full] <file> [file2]
```

(Equivalent explicit form, if you prefer: `java -cp <jar> uk.ac.ebi.ena.readtools.v2.cli.ValidateCli ...`)

### Options

| Option | Description |
|---|---|
| `--format`, `-f` | **Required.** `FASTQ`, `BAM`, or `CRAM`. |
| `--full` | Validate up to 100M reads (`webin-cli` "extended" mode). Default is **quick**: first 100k reads. |
| `--help`, `-h` | Show usage. |
| `<file> [file2]` | One file, or **two** FASTQ files for a paired submission. |

### Exit codes

| Code | Meaning |
|---|---|
| `0` | Valid |
| `1` | Invalid (validation error) |
| `2` | Usage error / file not found |

---

## Examples

Valid single-end FASTQ:

```
$ java -jar readtools-2.15.1-all.jar --format FASTQ reads.fastq.gz
RESULT: VALID
  format:       FASTQ
  mode:         quick (<=100k reads)
  reads.fastq.gz: reads=3, highQuality(avg>=30)=0
```

Valid paired-end FASTQ (pairing is detected and reported):

```
$ java -jar readtools-2.15.1-all.jar --format FASTQ R1.fastq.gz R2.fastq.gz
RESULT: VALID
  format:       FASTQ
  mode:         quick (<=100k reads)
  paired:       true
  R1.fastq.gz: reads=3829, highQuality(avg>=30)=3595
  R2.fastq.gz: reads=3829, highQuality(avg>=30)=2908
```

Invalid inputs (exit code 1):

```
$ ... --format FASTQ actually-a-fasta.fastq.gz
RESULT: INVALID
  Sequence header must start with @: >seq1 at line 1 in fastq

$ ... --format FASTQ non-iupac-base.fastq.gz
RESULT: INVALID
  Reads must contain only valid IUPAC codes

$ ... --format FASTQ length-mismatch.fastq.gz
RESULT: INVALID
  Sequence and quality line must be the same length at line 1 in fastq
```

---

## What it checks

The same rules readtools applies for ENA read submission, including:

- **Structure** (via htsjdk): each record starts with `@`, `+` separator line, bases/qualities
  present and equal length.
- **Bases**: valid IUPAC codes (`ACGTURYSWKMBDHVN.-`); rejects files where >50% of bases are
  non-AUTCG.
- **Quality**: printable ASCII; Phred offset (33/64) auto-detected and normalised.
- **Read names**: ≤256 chars; no duplicate read names; Casava 1.8 or `name/1`-style handling.
- **Paired FASTQ** (two files): pairing detected via read-name matching; fails if <20% of reads
  pair, or if naming doesn't allow pairing to be determined.

Compressed input (`.gz`, `.bz2`) and uncompressed FASTQ are all read transparently.

Note: only the first **100k** reads (quick) or **100M** reads (`--full`) are inspected; larger
files are validated partially.

---

## How the executable jar is wired

The fat jar's `Main-Class` is set in `build.gradle` so `java -jar` works directly:

```gradle
shadowJar {
    manifest {
        attributes( 'Main-Class': 'uk.ac.ebi.ena.readtools.v2.cli.ValidateCli' )
    }
}
```

For a cleaner upstream contribution you could instead wire `ValidateCli` in as a `validate`
subcommand of readtools' existing jcommander CLI in `cram/CramTools.java`, matching the
established pattern.
