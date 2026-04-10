# cutadapt-rs

A Rust reimplementation of [cutadapt](https://github.com/marcelm/cutadapt), the widely-used adapter trimming tool for high-throughput sequencing reads.

## Overview

cutadapt-rs is a 1:1 port of the Python cutadapt tool, providing the same functionality with the performance benefits of Rust. It removes adapter sequences, primers, poly-A tails, and other unwanted sequence from sequencing reads.

### Features

- **Adapter trimming** — Remove 3', 5', anchored, linked, and anywhere adapters
- **Quality trimming** — BWA-style quality trimming, NextSeq polyG trimming
- **Poly-A/T trimming** — Automatic poly-A tail removal
- **Read filtering** — Filter by length, N content, expected errors, CASAVA status
- **IUPAC wildcards** — Full IUPAC ambiguity code support in adapters and reads
- **Multiple adapters** — Specify multiple adapters; best match is selected
- **Paired-end support** — Process paired-end reads with coordinated filtering
- **Demultiplexing** — Route reads to different output files by adapter name
- **Multiple output actions** — Trim, mask with N, lowercase, retain adapter, crop to adapter

## Project Structure

```
cutadapt-rs/
  Cargo.toml                    # Workspace root
  crates/
    cutadapt-core/              # Library crate — all algorithms
      src/
        align.rs                # Semi-global DP alignment engine
        kmer_finder.rs          # Shift-and bitmatrix k-mer pre-screening
        kmer_heuristic.rs       # K-mer set generation for pre-screening
        qualtrim.rs             # Quality trimming (BWA, NextSeq, poly-A)
        match_tables.rs         # IUPAC character encoding tables
        record.rs               # SequenceRecord type
        info.rs                 # Per-read pipeline state
        adapters/               # Adapter matching hierarchy
          mod.rs                # Where flags, MultipleAdapters
          single.rs             # 9 adapter types (Front/Back/Anywhere/Prefix/Suffix/...)
          matches.rs            # RemoveBeforeMatch, RemoveAfterMatch, LinkedMatch
          linked.rs             # LinkedAdapter (front + back combination)
          statistics.rs         # Per-adapter match statistics
        modifiers/              # Read modification operations
          mod.rs                # SingleEndModifier / PairedEndModifier traits
          quality.rs            # QualityTrimmer, NextseqQualityTrimmer, PolyATrimmer
          simple.rs             # UnconditionalCutter, NEndTrimmer, Shortener, ZeroCapper
          rename.rs             # Renamer, PrefixSuffixAdder, LengthTagModifier
        pipeline.rs             # SingleEndPipeline, PairedEndPipeline
        steps.rs                # Filters, sinks, demultiplexers
        predicates.rs           # TooShort, TooLong, TooManyN, etc.
        parser.rs               # Adapter specification parsing
        files.rs                # FASTQ/FASTA I/O, compression detection
        report.rs               # Statistics collection and reporting
        json_output.rs          # JSON formatting helpers
        statistics.rs           # ReadLengthStatistics
        tokenizer.rs            # Template tokenization for --rename
    cutadapt/                   # Binary crate — CLI
      src/
        main.rs                 # clap argument parser (~100 options)
```

## Building

```bash
cargo build --release
```

## Usage

```bash
# Basic 3' adapter trimming
cutadapt -a AGATCGGAAGAGC -o output.fastq input.fastq

# 5' and 3' adapter trimming with quality filtering
cutadapt -g AATGATACGGCGACCACCGA -a AGATCGGAAGAGC -q 20 -m 25 -o output.fastq input.fastq

# Paired-end trimming
cutadapt -a ADAPT1 -A ADAPT2 -o out.1.fastq -p out.2.fastq in.1.fastq in.2.fastq

# Poly-A trimming
cutadapt --poly-a -o output.fastq input.fastq

# Multiple adapters, mask instead of trim
cutadapt -a ADAPTER1 -a ADAPTER2 --action=mask -o output.fastq input.fastq
```

Run `cutadapt --help` for the full list of options.

## Algorithm

The core alignment algorithm is a semi-global dynamic programming approach with:

- **Cost + score hybrid tracking** — Edit distance for error counting, alignment score for match quality
- **Ukkonen's trick** ��� Dynamic column pruning for O(kn) average-case complexity
- **K-mer pre-screening** — Shift-and bitmatrix algorithm filters ~99% of reads without running the full aligner
- **IUPAC wildcard support** — Bitwise AND matching for ambiguity codes

## Benchmark

Tested on 9,880,559 paired-end reads (150 bp, Illumina), standard TruSeq adapter trimming.

### Baseline

Commit: `0ad2e08`

| | Python cutadapt 5.2 | cutadapt-rs |
|---|---|---|
| Time | 1m 07s | 1m 47s |
| Reads processed | 9,880,559 | 9,880,559 |
| R1 with adapter | 306,997 (3.1%) | 306,997 (3.1%) |
| R2 with adapter | 298,514 (3.0%) | 298,514 (3.0%) |
| Total bp written | 2,961,299,044 | 2,961,299,044 |

### Current main

Measured with `./test.sh` on the same workload after the streaming pipeline, batched output, hot-path allocation cleanup, multi-threaded read processing, parallel gzip output, and `flate2` backend tuning.

#### thread=1

| Metric | Python cutadapt 5.2 | cutadapt-rs |
|---|---|---|
| Time | 1:02.55 | 0:48.33 |
| Reads processed | 9,880,559 | 9,880,559 |
| R1 with adapter | 306,997 (3.1%) | 306,997 (3.1%) |
| R2 with adapter | 298,514 (3.0%) | 298,514 (3.0%) |
| Total bp written | 2,961,299,044 bp | 2,961,299,044 bp |
| Peak RSS | n/a | ~60 MB |

#### thread=6

| Metric | Python cutadapt 5.2 | cutadapt-rs |
|---|---|---|
| Time | 0:19.33 | 0:21.11 |
| Reads processed | 9,880,559 | 9,880,559 |
| R1 with adapter | 306,997 (3.1%) | 306,997 (3.1%) |
| R2 with adapter | 298,514 (3.0%) | 298,514 (3.0%) |
| Total bp written | 2,961,299,044 bp | 2,961,299,044 bp |
| Peak RSS | n/a | ~226 MB |

Results remain numerically identical. Compared with the `0ad2e08` baseline, the main gains came from replacing `read_all()` with batched streaming, batching output writes, reusing modifiers outside the per-read loop, reducing hot-path cloning, parallelizing gzip output, parallelizing read processing, and tuning the gzip backend.

## Testing

```bash
# Run all tests (unit + integration)
cargo test

# Run only integration tests (golden file comparison with Python cutadapt)
cargo test --test integration_tests
```

The integration test suite compares Rust output byte-for-byte with expected outputs from Python cutadapt, covering adapter trimming, quality trimming, filtering, wildcards, multiple actions, and more.

## Compatibility

cutadapt-rs aims for exact behavioral compatibility with Python cutadapt. The same command-line arguments produce the same output. The integration test suite verifies this by comparing outputs against golden files from the Python test suite.

## License

Same license as the original cutadapt project (MIT).

## Acknowledgments

This is a Rust port of [cutadapt](https://github.com/marcelm/cutadapt) by Marcel Martin. All credit for the algorithms and design goes to the original author.
