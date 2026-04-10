use anyhow::{bail, Result};
use clap::Parser;
use std::fmt;

use cutadapt_core::files::{FileFormat, FileOpener, InputPaths};
use cutadapt_core::modifiers::adapter::AdapterTrimmer;
use cutadapt_core::modifiers::quality::{QualityTrimmer, NextseqQualityTrimmer};
use cutadapt_core::modifiers::simple::{UnconditionalCutter, NEndTrimmer, Shortener};
use cutadapt_core::modifiers::SingleEndModifier;
use cutadapt_core::parser::{make_adapters_from_specifications, AdapterType, SearchParameters};
use cutadapt_core::steps::PairFilterMode;
// ---------------------------------------------------------------------------
// Helper types
// ---------------------------------------------------------------------------

/// The action to take when an adapter match is found.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum MatchAction {
    Trim,
    Retain,
    Mask,
    Lowercase,
    Crop,
    None,
}

impl fmt::Display for MatchAction {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Trim => write!(f, "trim"),
            Self::Retain => write!(f, "retain"),
            Self::Mask => write!(f, "mask"),
            Self::Lowercase => write!(f, "lowercase"),
            Self::Crop => write!(f, "crop"),
            Self::None => write!(f, "none"),
        }
    }
}

impl std::str::FromStr for MatchAction {
    type Err = String;
    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        match s {
            "trim" => Ok(Self::Trim),
            "retain" => Ok(Self::Retain),
            "mask" => Ok(Self::Mask),
            "lowercase" => Ok(Self::Lowercase),
            "crop" => Ok(Self::Crop),
            "none" => Ok(Self::None),
            other => Err(format!(
                "unknown action '{other}': expected trim, retain, mask, lowercase, crop, or none"
            )),
        }
    }
}

/// Report style.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ReportStyle {
    Full,
    Minimal,
}

impl fmt::Display for ReportStyle {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Full => write!(f, "full"),
            Self::Minimal => write!(f, "minimal"),
        }
    }
}

impl std::str::FromStr for ReportStyle {
    type Err = String;
    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        match s {
            "full" => Ok(Self::Full),
            "minimal" => Ok(Self::Minimal),
            other => Err(format!(
                "unknown report style '{other}': expected full or minimal"
            )),
        }
    }
}

/// Pair-filter mode for paired-end filtering.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum PairFilter {
    Any,
    Both,
    First,
}

impl fmt::Display for PairFilter {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Any => write!(f, "any"),
            Self::Both => write!(f, "both"),
            Self::First => write!(f, "first"),
        }
    }
}

impl std::str::FromStr for PairFilter {
    type Err = String;
    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        match s {
            "any" => Ok(Self::Any),
            "both" => Ok(Self::Both),
            "first" => Ok(Self::First),
            other => Err(format!(
                "unknown pair filter '{other}': expected any, both, or first"
            )),
        }
    }
}

/// Where (which end) the adapter is expected.
#[derive(Debug, Clone)]
enum AdapterPosition {
    Back,
    Front,
    Anywhere,
}

/// A single adapter specification together with its expected position.
#[derive(Debug, Clone)]
struct AdapterSpec {
    position: AdapterPosition,
    spec: String,
}

impl fmt::Display for AdapterSpec {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let label = match self.position {
            AdapterPosition::Back => "back",
            AdapterPosition::Front => "front",
            AdapterPosition::Anywhere => "anywhere",
        };
        write!(f, "{}:{}", label, self.spec)
    }
}

// ---------------------------------------------------------------------------
// CLI definition  (mirrors Python cutadapt get_argument_parser())
// ---------------------------------------------------------------------------

/// cutadapt-rs -- Remove adapter sequences from high-throughput sequencing reads.
///
/// Usage:
///     cutadapt -a ADAPTER [options] [-o output.fastq] input.fastq
///
/// For paired-end reads:
///     cutadapt -a ADAPT1 -A ADAPT2 [options] -o out1.fastq -p out2.fastq in1.fastq in2.fastq
///
/// Replace "ADAPTER" with the actual sequence of your 3' adapter. IUPAC wildcard
/// characters are supported. All reads from input.fastq will be written to
/// output.fastq with the adapter sequence removed. Adapter matching is
/// error-tolerant. Multiple adapter sequences can be given (use further -a
/// options), but only the best-matching adapter will be removed.
///
/// Input may also be in FASTA format. Compressed input and output is supported and
/// auto-detected from the file name (.gz, .xz, .bz2). Use the file name '-' for
/// standard input/output. Without the -o option, output is sent to standard output.
#[derive(Parser, Debug)]
#[command(
    name = "cutadapt",
    version,
    about = "Remove adapter sequences from high-throughput sequencing reads",
    after_help = "Citation:\n\n\
        Marcel Martin. Cutadapt removes adapter sequences from high-throughput\n\
        sequencing reads. EMBnet.Journal, 17(1):10-12, May 2011.\n\
        http://dx.doi.org/10.14806/ej.17.1.200\n\n\
        See https://cutadapt.readthedocs.io/ for full documentation."
)]
struct Cli {
    // ── System / general ────────────────────────────────────────────────
    /// Print debug log. Use twice to also print DP matrices.
    #[arg(long, action = clap::ArgAction::Count)]
    debug: u8,

    /// Number of CPU cores to use. Use 0 to auto-detect. Default: 1
    #[arg(short = 'j', long, default_value_t = 1)]
    cores: i32,

    /// GC content as a percentage (hidden)
    #[arg(long, default_value_t = 50.0, hide = true)]
    gc_content: f64,

    /// Buffer size for the reader process when running in parallel (hidden)
    #[arg(long, default_value_t = 4_000_000, hide = true)]
    buffer_size: usize,

    /// Enable adapter index creation (hidden; use --no-index to disable)
    #[arg(long = "no-index", default_value_t = true, action = clap::ArgAction::SetFalse, hide = true)]
    index: bool,

    // ── Finding adapters ────────────────────────────────────────────────
    /// Sequence of an adapter ligated to the 3' end (paired data: of the
    /// first read). The adapter and subsequent bases are trimmed. If a '$'
    /// character is appended ('anchoring'), the adapter is only found if
    /// it is a suffix of the read.
    #[arg(short = 'a', long = "adapter", action = clap::ArgAction::Append, value_name = "ADAPTER")]
    back_adapters: Vec<String>,

    /// Sequence of an adapter ligated to the 5' end (paired data: of the
    /// first read). The adapter and any preceding bases are trimmed.
    /// Partial matches at the 5' end are allowed. If a '^' character is
    /// prepended ('anchoring'), the adapter is only found if it is a
    /// prefix of the read.
    #[arg(short = 'g', long = "front", action = clap::ArgAction::Append, value_name = "ADAPTER")]
    front_adapters: Vec<String>,

    /// Sequence of an adapter that may be ligated to the 5' or 3' end
    /// (paired data: of the first read).
    #[arg(short = 'b', long = "anywhere", action = clap::ArgAction::Append, value_name = "ADAPTER")]
    anywhere_adapters: Vec<String>,

    /// Maximum allowed error rate (if 0 <= E < 1), or absolute number of
    /// errors for full-length adapter match (if E is an integer >= 1).
    /// Default: 0.1 (10%)
    #[arg(short = 'e', long = "error-rate", aliases = ["errors"], default_value_t = 0.1, value_name = "E")]
    error_rate: f64,

    /// Allow only mismatches in alignments (no indels).
    #[arg(long = "no-indels")]
    no_indels: bool,

    /// Remove up to COUNT adapters from each read. Default: 1
    #[arg(short = 'n', long, default_value_t = 1, value_name = "COUNT")]
    times: usize,

    /// Require MINLENGTH overlap between read and adapter for an adapter
    /// to be found. Default: 3
    #[arg(short = 'O', long, default_value_t = 3, value_name = "MINLENGTH")]
    overlap: usize,

    /// Interpret IUPAC wildcards in reads. Default: false
    #[arg(long = "match-read-wildcards", short = 'N')]
    match_read_wildcards: bool,

    /// Do not interpret IUPAC wildcards in adapters.
    #[arg(long = "no-match-adapter-wildcards")]
    no_match_adapter_wildcards: bool,

    /// What to do if a match was found.
    /// trim: trim adapter and up-/downstream sequence;
    /// retain: trim, but retain adapter;
    /// mask: replace with 'N' characters;
    /// lowercase: convert to lowercase;
    /// crop: trim up and downstream sequence;
    /// none: leave unchanged. Default: trim
    #[arg(long, default_value = "trim", value_name = "ACTION")]
    action: MatchAction,

    /// Check both the read and its reverse complement for adapter matches.
    /// If match is on reverse-complemented version, output that one.
    #[arg(long = "rc", alias = "revcomp")]
    reverse_complement: bool,

    // ── Additional read modifications ───────────────────────────────────
    /// Remove LEN bases from each read (or R1 if paired; use -U for R2).
    /// If LEN is positive, remove bases from the beginning.
    /// If LEN is negative, remove bases from the end.
    /// Can be used twice if LENs have different signs.
    /// Applied *before* adapter trimming.
    #[arg(short = 'u', long = "cut", action = clap::ArgAction::Append, value_name = "LEN", allow_hyphen_values = true)]
    cut: Vec<i32>,

    /// NextSeq-specific quality trimming (each read). Trims also dark
    /// cycles appearing as high-quality G bases.
    #[arg(long = "nextseq-trim", value_name = "3'CUTOFF")]
    nextseq_trim: Option<u8>,

    /// Trim low-quality bases from 5' and/or 3' ends of each read before
    /// adapter removal. Applied to both reads if data is paired. If one
    /// value is given, only the 3' end is trimmed. If two comma-separated
    /// cutoffs are given, the 5' end is trimmed with the first cutoff, the
    /// 3' end with the second.
    #[arg(short = 'q', long = "quality-cutoff", value_name = "[5'CUTOFF,]3'CUTOFF")]
    quality_cutoff: Option<String>,

    /// Assume that quality values in FASTQ are encoded as ascii(quality + N).
    /// This needs to be set to 64 for some old Illumina FASTQ files. Default: 33
    #[arg(long = "quality-base", default_value_t = 33, value_name = "N")]
    quality_base: u8,

    /// Trim poly-A tails.
    #[arg(long = "poly-a")]
    poly_a: bool,

    /// Shorten reads to LENGTH. Positive values remove bases at the end
    /// while negative ones remove bases at the beginning. Applied after
    /// adapter trimming.
    #[arg(short = 'l', long = "length", value_name = "LENGTH", allow_hyphen_values = true)]
    length: Option<i32>,

    /// Trim N's on ends of reads.
    #[arg(long = "trim-n")]
    trim_n: bool,

    /// Search for TAG followed by a decimal number in the description
    /// field of the read. Replace the decimal number with the correct
    /// length of the trimmed read.
    #[arg(long = "length-tag", value_name = "TAG")]
    length_tag: Option<String>,

    /// Remove this suffix from read names if present. Can be given
    /// multiple times.
    #[arg(long = "strip-suffix", action = clap::ArgAction::Append)]
    strip_suffix: Vec<String>,

    /// Add this prefix to read names. Use {name} to insert the name of the
    /// matching adapter.
    #[arg(short = 'x', long = "prefix", default_value = "")]
    prefix: String,

    /// Add this suffix to read names; can also include {name}.
    #[arg(short = 'y', long = "suffix", default_value = "")]
    suffix: String,

    /// Rename reads using TEMPLATE containing variables such as {id},
    /// {adapter_name} etc. (see documentation)
    #[arg(long, value_name = "TEMPLATE")]
    rename: Option<String>,

    /// Change negative quality values to zero.
    #[arg(short = 'z', long = "zero-cap")]
    zero_cap: bool,

    // ── Filtering of processed reads ────────────────────────────────────
    /// Discard reads shorter than LEN. Default: 0.
    /// For paired-end data, use LEN:LEN2 to set separate limits.
    #[arg(short = 'm', long = "minimum-length", value_name = "LEN[:LEN2]")]
    minimum_length: Option<String>,

    /// Discard reads longer than LEN. Default: no limit.
    /// For paired-end data, use LEN:LEN2 to set separate limits.
    #[arg(short = 'M', long = "maximum-length", value_name = "LEN[:LEN2]")]
    maximum_length: Option<String>,

    /// Discard reads with more than COUNT 'N' bases. If COUNT is a number
    /// between 0 and 1, it is interpreted as a fraction of the read length.
    #[arg(long = "max-n", value_name = "COUNT")]
    max_n: Option<f64>,

    /// Discard reads whose expected number of errors (computed from quality
    /// values) exceeds ERRORS.
    #[arg(long = "max-expected-errors", alias = "max-ee", value_name = "ERRORS")]
    max_expected_errors: Option<f64>,

    /// As --max-expected-errors, but divided by length to account for reads
    /// of varying length.
    #[arg(long = "max-average-error-rate", alias = "max-aer", value_name = "ERROR_RATE")]
    max_average_error_rate: Option<f64>,

    /// Discard reads that contain an adapter. Use also -O to avoid
    /// discarding too many randomly matching reads.
    #[arg(long = "discard-trimmed", alias = "discard")]
    discard_trimmed: bool,

    /// Discard reads that do not contain an adapter.
    #[arg(long = "discard-untrimmed", alias = "trimmed-only")]
    discard_untrimmed: bool,

    /// Discard reads that did not pass CASAVA filtering (header has :Y:).
    #[arg(long = "discard-casava")]
    discard_casava: bool,

    // ── Output ──────────────────────────────────────────────────────────
    /// Print only error messages.
    #[arg(long)]
    quiet: bool,

    /// Which type of report to print: 'full' or 'minimal'. Default: full
    #[arg(long, value_name = "STYLE")]
    report: Option<ReportStyle>,

    /// Dump report in JSON format to FILE.
    #[arg(long = "json", value_name = "FILE")]
    json: Option<String>,

    /// Write trimmed reads to FILE. FASTQ or FASTA format is chosen
    /// depending on input. Use '{name}' for demultiplexing (see docs).
    /// Default: write to standard output
    #[arg(short = 'o', long = "output", value_name = "FILE")]
    output: Option<String>,

    /// Output FASTA to standard output even on FASTQ input.
    #[arg(long)]
    fasta: bool,

    /// Compression level for compressed output files. Default: 1
    #[arg(long = "compression-level", default_value_t = 1, value_name = "N")]
    compression_level: u32,

    /// Write information about each read and its adapter matches into FILE.
    #[arg(long = "info-file", value_name = "FILE")]
    info_file: Option<String>,

    /// When the adapter matches in the middle of a read, write the rest
    /// (after the adapter) to FILE.
    #[arg(short = 'r', long = "rest-file", value_name = "FILE")]
    rest_file: Option<String>,

    /// When the adapter has N wildcard bases, write adapter bases matching
    /// wildcard positions to FILE.
    #[arg(long = "wildcard-file", value_name = "FILE")]
    wildcard_file: Option<String>,

    /// Write reads that are too short (according to -m) to FILE.
    #[arg(long = "too-short-output", value_name = "FILE")]
    too_short_output: Option<String>,

    /// Write reads that are too long (according to -M) to FILE.
    #[arg(long = "too-long-output", value_name = "FILE")]
    too_long_output: Option<String>,

    /// Write reads that do not contain any adapter to FILE.
    #[arg(long = "untrimmed-output", value_name = "FILE")]
    untrimmed_output: Option<String>,

    // ── Paired-end options ──────────────────────────────────────────────
    /// 3' adapter to be removed from R2.
    #[arg(short = 'A', action = clap::ArgAction::Append, value_name = "ADAPTER")]
    back_adapters2: Vec<String>,

    /// 5' adapter to be removed from R2.
    #[arg(short = 'G', action = clap::ArgAction::Append, value_name = "ADAPTER")]
    front_adapters2: Vec<String>,

    /// 5'/3' adapter to be removed from R2.
    #[arg(short = 'B', action = clap::ArgAction::Append, value_name = "ADAPTER")]
    anywhere_adapters2: Vec<String>,

    /// Remove LENGTH bases from R2.
    #[arg(short = 'U', action = clap::ArgAction::Append, value_name = "LENGTH", allow_hyphen_values = true)]
    cut2: Vec<i32>,

    /// Quality-trimming cutoff for R2. Default: same as for R1.
    #[arg(short = 'Q', long = "quality-cutoff2", value_name = "[5'CUTOFF,]3'CUTOFF")]
    quality_cutoff2: Option<String>,

    /// Shorten R2 to LENGTH. Default: same as for R1.
    #[arg(short = 'L', value_name = "LENGTH", allow_hyphen_values = true)]
    length2: Option<i32>,

    /// Write R2 to FILE.
    #[arg(short = 'p', long = "paired-output", value_name = "FILE")]
    paired_output: Option<String>,

    /// Write info about R2 to FILE (see --info-file).
    #[arg(long = "info-file-paired", value_name = "FILE")]
    info_file2: Option<String>,

    /// Treat adapters given with -a/-A etc. as pairs. Either both or none
    /// are removed from each read pair.
    #[arg(long = "pair-adapters")]
    pair_adapters: bool,

    /// Which of the reads in a paired-end read have to match the filtering
    /// criterion in order for the pair to be filtered. Default: any
    #[arg(long = "pair-filter", value_name = "MODE")]
    pair_filter: Option<PairFilter>,

    /// Read and/or write interleaved paired-end reads.
    #[arg(long)]
    interleaved: bool,

    /// Write second read in a pair to this FILE when no adapter was found.
    /// Use with --untrimmed-output.
    #[arg(long = "untrimmed-paired-output", value_name = "FILE")]
    untrimmed_paired_output: Option<String>,

    /// Write second read in a pair to this file if pair is too short.
    #[arg(long = "too-short-paired-output", value_name = "FILE")]
    too_short_paired_output: Option<String>,

    /// Write second read in a pair to this file if pair is too long.
    #[arg(long = "too-long-paired-output", value_name = "FILE")]
    too_long_paired_output: Option<String>,

    // ── Positional ──────────────────────────────────────────────────────
    /// Input FASTQ/FASTA file(s). Provide one for single-end, two for
    /// paired-end reads.
    #[arg(value_name = "INPUT")]
    inputs: Vec<String>,
}

// ---------------------------------------------------------------------------
// Helper functions mirroring the Python ones
// ---------------------------------------------------------------------------

/// Parse a quality-cutoff string of the form `INT` or `INT,INT` into
/// `(five_prime_cutoff, three_prime_cutoff)`.
///
/// A single value is treated as a 3'-only cutoff (5' cutoff = 0).
///
/// ```text
/// "5"   -> (0, 5)
/// "6,7" -> (6, 7)
/// ```
fn parse_cutoffs(s: &str) -> Result<(u8, u8)> {
    let parts: Vec<&str> = s.split(',').collect();
    match parts.len() {
        1 => {
            let v: u8 = parts[0]
                .parse()
                .map_err(|e| anyhow::anyhow!("Quality cutoff value not recognised: {e}"))?;
            Ok((0, v))
        }
        2 => {
            let a: u8 = parts[0]
                .parse()
                .map_err(|e| anyhow::anyhow!("Quality cutoff value not recognised: {e}"))?;
            let b: u8 = parts[1]
                .parse()
                .map_err(|e| anyhow::anyhow!("Quality cutoff value not recognised: {e}"))?;
            Ok((a, b))
        }
        _ => bail!(
            "Expected one value or two values separated by comma for the quality cutoff"
        ),
    }
}

/// Parse a length string of the form `INT`, `INT:INT`, `INT:`, or `:INT` into
/// a vector of `Option<i64>` values.
///
/// ```text
/// "25"    -> [Some(25)]
/// "17:25" -> [Some(17), Some(25)]
/// "25:"   -> [Some(25), None]
/// ":25"   -> [None,     Some(25)]
/// ```
fn parse_lengths(s: &str) -> Result<Vec<Option<i64>>> {
    let fields: Vec<&str> = s.split(':').collect();
    if fields.len() > 2 {
        bail!("Only at most one colon is allowed");
    }
    let values: Vec<Option<i64>> = fields
        .iter()
        .map(|f| {
            if f.is_empty() {
                Ok(None)
            } else {
                f.parse::<i64>()
                    .map(Some)
                    .map_err(|e| anyhow::anyhow!("Value not recognised: {e}"))
            }
        })
        .collect::<Result<Vec<_>>>()?;

    if values.len() == 2 && values[0].is_none() && values[1].is_none() {
        bail!("Cannot parse '{s}': At least one length needs to be given");
    }
    Ok(values)
}

// ---------------------------------------------------------------------------
// Collect adapter specs from the three flag lists
// ---------------------------------------------------------------------------

impl Cli {
    /// Merge `-a`, `-g`, `-b` into a single `Vec<AdapterSpec>` for R1.
    fn adapters(&self) -> Vec<AdapterSpec> {
        let mut out = Vec::new();
        for s in &self.back_adapters {
            out.push(AdapterSpec {
                position: AdapterPosition::Back,
                spec: s.clone(),
            });
        }
        for s in &self.front_adapters {
            out.push(AdapterSpec {
                position: AdapterPosition::Front,
                spec: s.clone(),
            });
        }
        for s in &self.anywhere_adapters {
            out.push(AdapterSpec {
                position: AdapterPosition::Anywhere,
                spec: s.clone(),
            });
        }
        out
    }

    /// Merge `-A`, `-G`, `-B` into a single `Vec<AdapterSpec>` for R2.
    fn adapters2(&self) -> Vec<AdapterSpec> {
        let mut out = Vec::new();
        for s in &self.back_adapters2 {
            out.push(AdapterSpec {
                position: AdapterPosition::Back,
                spec: s.clone(),
            });
        }
        for s in &self.front_adapters2 {
            out.push(AdapterSpec {
                position: AdapterPosition::Front,
                spec: s.clone(),
            });
        }
        for s in &self.anywhere_adapters2 {
            out.push(AdapterSpec {
                position: AdapterPosition::Anywhere,
                spec: s.clone(),
            });
        }
        out
    }

    /// Whether we should operate in paired-end mode (mirrors Python
    /// `determine_paired`).
    fn is_paired(&self) -> bool {
        self.paired_output.is_some()
            || self.interleaved
            || !self.adapters2().is_empty()
            || !self.cut2.is_empty()
            || self.length2.is_some()
            || self.pair_filter.is_some()
            || self.untrimmed_paired_output.is_some()
            || self.too_short_paired_output.is_some()
            || self.too_long_paired_output.is_some()
            || self.quality_cutoff2.is_some()
            || self.info_file2.is_some()
    }
}

// ---------------------------------------------------------------------------
// Pretty-print parsed configuration
// ---------------------------------------------------------------------------

fn print_config(cli: &Cli) {
    let paired = cli.is_paired();
    let adapters = cli.adapters();
    let adapters2 = cli.adapters2();

    println!("cutadapt-rs (Rust port)");
    println!("=======================");
    println!();

    // -- Input --
    println!("Input files:");
    if cli.inputs.is_empty() {
        println!("  (none — will read from stdin)");
    } else {
        for (i, path) in cli.inputs.iter().enumerate() {
            println!("  [{}] {}", i + 1, path);
        }
    }
    println!("  Paired-end mode: {}", paired);
    println!("  Interleaved:     {}", cli.interleaved);
    println!();

    // -- Adapters R1 --
    if !adapters.is_empty() {
        println!("R1 adapters ({}):", adapters.len());
        for a in &adapters {
            println!("  {}", a);
        }
        println!();
    }

    // -- Adapters R2 --
    if !adapters2.is_empty() {
        println!("R2 adapters ({}):", adapters2.len());
        for a in &adapters2 {
            println!("  {}", a);
        }
        println!();
    }

    // -- Adapter parameters --
    println!("Adapter matching parameters:");
    println!("  Error rate:              {}", cli.error_rate);
    println!("  Times:                   {}", cli.times);
    println!("  Overlap:                 {}", cli.overlap);
    println!("  Match read wildcards:    {}", cli.match_read_wildcards);
    println!("  Match adapter wildcards: {}", !cli.no_match_adapter_wildcards);
    println!("  Indels allowed:          {}", !cli.no_indels);
    println!("  Action:                  {}", cli.action);
    println!("  Reverse complement:      {}", cli.reverse_complement);
    if cli.pair_adapters {
        println!("  Pair adapters:           true");
    }
    println!();

    // -- Additional modifications --
    println!("Additional modifications:");
    if !cli.cut.is_empty() {
        println!("  Unconditional cut (R1):  {:?}", cli.cut);
    }
    if !cli.cut2.is_empty() {
        println!("  Unconditional cut (R2):  {:?}", cli.cut2);
    }
    if let Some(v) = cli.nextseq_trim {
        println!("  NextSeq trim:            {}", v);
    }
    if let Some(ref q) = cli.quality_cutoff {
        match parse_cutoffs(q) {
            Ok((five, three)) => println!("  Quality cutoff (R1):     5'={}, 3'={}", five, three),
            Err(e) => println!("  Quality cutoff (R1):     ERROR: {}", e),
        }
    }
    if let Some(ref q2) = cli.quality_cutoff2 {
        match parse_cutoffs(q2) {
            Ok((five, three)) => println!("  Quality cutoff (R2):     5'={}, 3'={}", five, three),
            Err(e) => println!("  Quality cutoff (R2):     ERROR: {}", e),
        }
    }
    println!("  Quality base:            {}", cli.quality_base);
    println!("  Poly-A trimming:         {}", cli.poly_a);
    if let Some(len) = cli.length {
        println!("  Shorten reads (R1):      {}", len);
    }
    if let Some(len) = cli.length2 {
        println!("  Shorten reads (R2):      {}", len);
    }
    println!("  Trim N:                  {}", cli.trim_n);
    if let Some(ref tag) = cli.length_tag {
        println!("  Length tag:              {}", tag);
    }
    if !cli.strip_suffix.is_empty() {
        println!("  Strip suffixes:          {:?}", cli.strip_suffix);
    }
    if !cli.prefix.is_empty() {
        println!("  Prefix:                  {}", cli.prefix);
    }
    if !cli.suffix.is_empty() {
        println!("  Suffix:                  {}", cli.suffix);
    }
    if let Some(ref tpl) = cli.rename {
        println!("  Rename template:         {}", tpl);
    }
    println!("  Zero cap:                {}", cli.zero_cap);
    println!();

    // -- Filtering --
    println!("Filtering:");
    if let Some(ref ml) = cli.minimum_length {
        match parse_lengths(ml) {
            Ok(vals) => println!("  Minimum length:          {:?}", vals),
            Err(e) => println!("  Minimum length:          ERROR: {}", e),
        }
    } else {
        println!("  Minimum length:          0 (default)");
    }
    if let Some(ref ml) = cli.maximum_length {
        match parse_lengths(ml) {
            Ok(vals) => println!("  Maximum length:          {:?}", vals),
            Err(e) => println!("  Maximum length:          ERROR: {}", e),
        }
    } else {
        println!("  Maximum length:          (no limit)");
    }
    if let Some(v) = cli.max_n {
        println!("  Max N:                   {}", v);
    }
    if let Some(v) = cli.max_expected_errors {
        println!("  Max expected errors:     {}", v);
    }
    if let Some(v) = cli.max_average_error_rate {
        println!("  Max average error rate:  {}", v);
    }
    println!("  Discard trimmed:         {}", cli.discard_trimmed);
    println!("  Discard untrimmed:       {}", cli.discard_untrimmed);
    println!("  Discard CASAVA:          {}", cli.discard_casava);
    if let Some(pf) = cli.pair_filter {
        println!("  Pair filter:             {}", pf);
    }
    println!();

    // -- Output --
    println!("Output:");
    println!("  Quiet:                   {}", cli.quiet);
    if let Some(rs) = cli.report {
        println!("  Report:                  {}", rs);
    } else {
        println!("  Report:                  full (default)");
    }
    if let Some(ref p) = cli.json {
        println!("  JSON report:             {}", p);
    }
    if let Some(ref p) = cli.output {
        println!("  Output (R1):             {}", p);
    } else {
        println!("  Output (R1):             (stdout)");
    }
    if let Some(ref p) = cli.paired_output {
        println!("  Output (R2):             {}", p);
    }
    println!("  Force FASTA:             {}", cli.fasta);
    println!("  Compression level:       {}", cli.compression_level);
    if let Some(ref p) = cli.info_file {
        println!("  Info file:               {}", p);
    }
    if let Some(ref p) = cli.info_file2 {
        println!("  Info file (R2):          {}", p);
    }
    if let Some(ref p) = cli.rest_file {
        println!("  Rest file:               {}", p);
    }
    if let Some(ref p) = cli.wildcard_file {
        println!("  Wildcard file:           {}", p);
    }
    if let Some(ref p) = cli.too_short_output {
        println!("  Too-short output:        {}", p);
    }
    if let Some(ref p) = cli.too_short_paired_output {
        println!("  Too-short paired output: {}", p);
    }
    if let Some(ref p) = cli.too_long_output {
        println!("  Too-long output:         {}", p);
    }
    if let Some(ref p) = cli.too_long_paired_output {
        println!("  Too-long paired output:  {}", p);
    }
    if let Some(ref p) = cli.untrimmed_output {
        println!("  Untrimmed output:        {}", p);
    }
    if let Some(ref p) = cli.untrimmed_paired_output {
        println!("  Untrimmed paired output: {}", p);
    }
    println!();

    // -- System --
    println!("System:");
    println!("  Cores:                   {}", cli.cores);
    println!("  Buffer size:             {}", cli.buffer_size);
    println!("  Debug level:             {}", cli.debug);
    println!("  Index enabled:           {}", cli.index);
}

// ---------------------------------------------------------------------------
// Entry point
// ---------------------------------------------------------------------------

fn main() -> Result<()> {
    let cli = Cli::parse();

    // Initialise logging.
    let log_level = match cli.debug {
        0 => {
            if cli.quiet {
                log::LevelFilter::Error
            } else {
                log::LevelFilter::Info
            }
        }
        1 => log::LevelFilter::Debug,
        _ => log::LevelFilter::Trace,
    };
    env_logger::Builder::new().filter_level(log_level).init();

    // --- Validate a few things early (mirrors check_arguments) ----------
    if cli.overlap < 1 {
        bail!("The overlap must be at least 1.");
    }
    if cli.cores < 0 {
        bail!("Value for --cores cannot be negative");
    }
    if !(0.0..=100.0).contains(&cli.gc_content) {
        bail!("GC content must be given as percentage between 0 and 100");
    }
    if cli.pair_adapters && cli.times != 1 {
        bail!("--pair-adapters cannot be used with --times");
    }
    if cli.quiet && cli.report.is_some() {
        bail!("Options --quiet and --report cannot be used at the same time");
    }
    if cli.rename.is_some() && (!cli.prefix.is_empty() || !cli.suffix.is_empty()) {
        bail!("Option --rename cannot be combined with --prefix (-x) or --suffix (-y)");
    }

    let paired = cli.is_paired();
    if paired && !cli.interleaved {
        if cli.paired_output.is_none() {
            bail!(
                "When a paired-end trimming option such as -A/-G/-B/-U is used, \
                 a second output file needs to be specified via -p (--paired-output)."
            );
        }
        if cli.output.is_none() {
            bail!(
                "When you use -p or --paired-output, you must also use the -o option."
            );
        }
    }

    // Show what we parsed.
<<<<<<< HEAD
    print_config(&cli);

    println!("NOTE: Actual read processing is not yet implemented.");
    println!("      The pipeline will be wired up in a future step.");
=======
    if !cli.quiet {
        print_config(&cli);
    }

    run_pipeline(&cli)
}

fn run_pipeline(cli: &Cli) -> Result<()> {
    let paired = cli.is_paired();

    // ── Build search parameters ──────────────────────────────────────────
    let search_params = SearchParameters {
        max_errors: cli.error_rate,
        min_overlap: cli.overlap,
        read_wildcards: cli.match_read_wildcards,
        adapter_wildcards: !cli.no_match_adapter_wildcards,
        indels: !cli.no_indels,
    };

    // ── Build adapters ───────────────────────────────────────────────────
    let adapter_specs1: Vec<(AdapterType, &str)> = cli
        .back_adapters.iter().map(|s| (AdapterType::Back, s.as_str()))
        .chain(cli.front_adapters.iter().map(|s| (AdapterType::Front, s.as_str())))
        .chain(cli.anywhere_adapters.iter().map(|s| (AdapterType::Anywhere, s.as_str())))
        .collect();

    let adapter_specs2: Vec<(AdapterType, &str)> = cli
        .back_adapters2.iter().map(|s| (AdapterType::Back, s.as_str()))
        .chain(cli.front_adapters2.iter().map(|s| (AdapterType::Front, s.as_str())))
        .chain(cli.anywhere_adapters2.iter().map(|s| (AdapterType::Anywhere, s.as_str())))
        .collect();

    let adapters1 = make_adapters_from_specifications(&adapter_specs1, &search_params)
        .map_err(|e| anyhow::anyhow!("Adapter parse error (R1): {e}"))?;
    let adapters2 = make_adapters_from_specifications(&adapter_specs2, &search_params)
        .map_err(|e| anyhow::anyhow!("Adapter parse error (R2): {e}"))?;

    // ── Read input ───────────────────────────────────────────────────────
    let opener = FileOpener::new(cli.compression_level as i32, None);

    let input_paths = if paired && !cli.interleaved && cli.inputs.len() == 2 {
        InputPaths::new(cli.inputs.clone(), false)
    } else {
        InputPaths::new(cli.inputs.clone(), cli.interleaved)
    };

    let mut input_files = input_paths.open()?;
    let format = input_files.format;
    let has_qualities = format.has_qualities();

    let records1 = input_files.read_all()?;
    let records2 = if paired && !cli.interleaved && input_files.readers.len() > 1 {
        let reader = &mut input_files.readers[1];
        match format {
            FileFormat::Fastq => cutadapt_core::files::read_fastq(reader)?,
            FileFormat::Fasta => cutadapt_core::files::read_fasta(reader)?,
        }
    } else {
        Vec::new()
    };

    // ── Quality cutoffs ──────────────────────────────────────────────────
    let quality_base = cli.quality_base as i32;

    let qual_cutoff1 = cli.quality_cutoff.as_deref()
        .map(parse_cutoffs)
        .transpose()?;
    let qual_cutoff2 = cli.quality_cutoff2.as_deref()
        .map(parse_cutoffs)
        .transpose()?
        .or(qual_cutoff1);

    // ── Pair filter mode ─────────────────────────────────────────────────
    let pair_filter_mode = match cli.pair_filter {
        Some(PairFilter::Both) => PairFilterMode::Both,
        Some(PairFilter::First) => PairFilterMode::First,
        _ => PairFilterMode::Any,
    };

    // ── Min/max lengths ──────────────────────────────────────────────────
    let (min_len1, min_len2) = parse_paired_length(cli.minimum_length.as_deref())?;
    let (max_len1, max_len2) = parse_paired_length(cli.maximum_length.as_deref())?;

    // ── Run pipeline ─────────────────────────────────────────────────────
    let (n_reads, n_written, bp_in1, bp_in2, bp_out1, bp_out2,
         reads_with_adapter1, reads_with_adapter2) = if paired {
        run_paired(
            cli, records1, records2, adapters1, adapters2,
            &search_params, &opener, format,
            qual_cutoff1, qual_cutoff2, quality_base,
            pair_filter_mode, min_len1, min_len2, max_len1, max_len2,
        )?
    } else {
        run_single_direct(
            cli, &opener, format,
            qual_cutoff1, quality_base,
            min_len1, max_len1,
        )?
    };

    // ── Print report ─────────────────────────────────────────────────────
    if !cli.quiet {
        print_report(
            paired, n_reads, n_written, bp_in1, bp_in2, bp_out1, bp_out2,
            reads_with_adapter1, reads_with_adapter2,
        );
    }
>>>>>>> 0ad2e08 (Wire up read processing pipeline in main.rs)

    Ok(())
}

<<<<<<< HEAD
=======
#[allow(clippy::too_many_arguments)]
fn run_single_direct(
    cli: &Cli,
    opener: &FileOpener,
    format: FileFormat,
    qual_cutoff: Option<(u8, u8)>,
    quality_base: i32,
    min_len: Option<usize>,
    max_len: Option<usize>,
) -> Result<(u64, u64, u64, u64, u64, u64, u64, u64)> {
    let input_paths = InputPaths::new(cli.inputs.clone(), cli.interleaved);
    let mut input_files = input_paths.open()?;
    let records = input_files.read_all()?;
    let bp_in1: u64 = records.iter().map(|r| r.len() as u64).sum();
    let n_reads = records.len() as u64;

    let search_params = SearchParameters {
        max_errors: cli.error_rate,
        min_overlap: cli.overlap,
        read_wildcards: cli.match_read_wildcards,
        adapter_wildcards: !cli.no_match_adapter_wildcards,
        indels: !cli.no_indels,
    };
    let adapter_specs1: Vec<(AdapterType, &str)> = cli
        .back_adapters.iter().map(|s| (AdapterType::Back, s.as_str()))
        .chain(cli.front_adapters.iter().map(|s| (AdapterType::Front, s.as_str())))
        .chain(cli.anywhere_adapters.iter().map(|s| (AdapterType::Anywhere, s.as_str())))
        .collect();
    let adapters1 = make_adapters_from_specifications(&adapter_specs1, &search_params)
        .map_err(|e| anyhow::anyhow!("{e}"))?;

    let mut trimmer = AdapterTrimmer::new(adapters1, cli.times);

    let output_path = cli.output.as_deref().unwrap_or("-");
    let mut writer = opener.open_write(output_path)?;

    let mut n_written: u64 = 0;
    let mut bp_out1: u64 = 0;

    for record in records {
        let mut info = cutadapt_core::info::ModificationInfo::new(record.clone());
        let mut r = record;

        // Pre-adapter cuts
        for &cut in &cli.cut {
            let mut cutter = UnconditionalCutter::new(cut);
            r = cutter.modify(r, &mut info);
        }
        // Quality trimming
        if let Some(cutoff) = cli.nextseq_trim {
            let mut qt = NextseqQualityTrimmer::new(cutoff as i32, quality_base);
            r = qt.modify(r, &mut info);
        } else if let Some((five, three)) = qual_cutoff {
            if five > 0 || three > 0 {
                let mut qt = QualityTrimmer::new(five as i32, three as i32, quality_base);
                r = qt.modify(r, &mut info);
            }
        }
        // Adapter trimming
        r = trimmer.modify(r, &mut info);
        // Post-adapter
        if cli.trim_n {
            let mut nt = NEndTrimmer::new();
            r = nt.modify(r, &mut info);
        }
        if let Some(len) = cli.length {
            let mut sh = Shortener::new(len);
            r = sh.modify(r, &mut info);
        }

        // Filters
        if let Some(min) = min_len {
            if r.len() < min { continue; }
        }
        if let Some(max) = max_len {
            if r.len() > max { continue; }
        }
        if let Some(n) = cli.max_n {
            let n_count = r.sequence.chars().filter(|&c| c == 'N' || c == 'n').count();
            let threshold = if n < 1.0 { (n * r.len() as f64) as usize } else { n as usize };
            if n_count > threshold { continue; }
        }

        // Write
        bp_out1 += r.len() as u64;
        n_written += 1;
        write_record(&mut *writer, &r, format)?;
    }
    writer.flush()?;

    Ok((n_reads, n_written, bp_in1, 0, bp_out1, 0,
        trimmer.reads_with_adapter, 0))
}

#[allow(clippy::too_many_arguments)]
fn run_paired(
    cli: &Cli,
    records1: Vec<cutadapt_core::record::SequenceRecord>,
    records2: Vec<cutadapt_core::record::SequenceRecord>,
    adapters1: Vec<cutadapt_core::adapters::Adapter>,
    adapters2: Vec<cutadapt_core::adapters::Adapter>,
    _search_params: &SearchParameters,
    opener: &FileOpener,
    format: FileFormat,
    qual_cutoff1: Option<(u8, u8)>,
    qual_cutoff2: Option<(u8, u8)>,
    quality_base: i32,
    pair_filter_mode: PairFilterMode,
    min_len1: Option<usize>,
    min_len2: Option<usize>,
    max_len1: Option<usize>,
    max_len2: Option<usize>,
) -> Result<(u64, u64, u64, u64, u64, u64, u64, u64)> {
    let n_reads = records1.len() as u64;
    let bp_in1: u64 = records1.iter().map(|r| r.len() as u64).sum();
    let bp_in2: u64 = records2.iter().map(|r| r.len() as u64).sum();

    let mut trimmer1 = AdapterTrimmer::new(adapters1, cli.times);
    let mut trimmer2 = AdapterTrimmer::new(adapters2, cli.times);

    let out1_path = cli.output.as_deref().unwrap_or("-");
    let out2_path = cli.paired_output.as_deref().unwrap_or("-");
    let mut writer1 = opener.open_write(out1_path)?;
    let mut writer2 = opener.open_write(out2_path)?;

    let mut n_written: u64 = 0;
    let mut bp_out1: u64 = 0;
    let mut bp_out2: u64 = 0;

    for (r1_orig, r2_orig) in records1.into_iter().zip(records2.into_iter()) {
        let mut info1 = cutadapt_core::info::ModificationInfo::new(r1_orig.clone());
        let mut info2 = cutadapt_core::info::ModificationInfo::new(r2_orig.clone());
        let mut r1 = r1_orig;
        let mut r2 = r2_orig;

        // Pre-adapter cuts R1
        for &cut in &cli.cut {
            let mut c = UnconditionalCutter::new(cut);
            r1 = c.modify(r1, &mut info1);
        }
        // Pre-adapter cuts R2
        for &cut in &cli.cut2 {
            let mut c = UnconditionalCutter::new(cut);
            r2 = c.modify(r2, &mut info2);
        }
        // Quality trimming R1
        if let Some(cutoff) = cli.nextseq_trim {
            let mut qt = NextseqQualityTrimmer::new(cutoff as i32, quality_base);
            r1 = qt.modify(r1, &mut info1);
            r2 = qt.modify(r2, &mut info2);
        } else {
            if let Some((five, three)) = qual_cutoff1 {
                if five > 0 || three > 0 {
                    let mut qt = QualityTrimmer::new(five as i32, three as i32, quality_base);
                    r1 = qt.modify(r1, &mut info1);
                }
            }
            if let Some((five, three)) = qual_cutoff2 {
                if five > 0 || three > 0 {
                    let mut qt = QualityTrimmer::new(five as i32, three as i32, quality_base);
                    r2 = qt.modify(r2, &mut info2);
                }
            }
        }
        // Adapter trimming
        r1 = trimmer1.modify(r1, &mut info1);
        r2 = trimmer2.modify(r2, &mut info2);
        // Post-adapter
        if cli.trim_n {
            let mut nt = NEndTrimmer::new();
            r1 = nt.modify(r1, &mut info1);
            r2 = nt.modify(r2, &mut info2);
        }
        if let Some(len) = cli.length {
            let mut sh = Shortener::new(len);
            r1 = sh.modify(r1, &mut info1);
        }
        if let Some(len) = cli.length2.or(cli.length) {
            let mut sh = Shortener::new(len);
            r2 = sh.modify(r2, &mut info2);
        }

        // Paired filters
        let filter1 = check_filters(&r1, min_len1, max_len1, cli.max_n);
        let filter2 = check_filters(&r2, min_len2, max_len2, cli.max_n);
        let discard = match pair_filter_mode {
            PairFilterMode::Any => filter1 || filter2,
            PairFilterMode::Both => filter1 && filter2,
            PairFilterMode::First => filter1,
        };
        if discard { continue; }

        bp_out1 += r1.len() as u64;
        bp_out2 += r2.len() as u64;
        n_written += 1;
        write_record(&mut *writer1, &r1, format)?;
        write_record(&mut *writer2, &r2, format)?;
    }
    writer1.flush()?;
    writer2.flush()?;

    Ok((n_reads, n_written, bp_in1, bp_in2, bp_out1, bp_out2,
        trimmer1.reads_with_adapter, trimmer2.reads_with_adapter))
}

fn check_filters(
    r: &cutadapt_core::record::SequenceRecord,
    min_len: Option<usize>,
    max_len: Option<usize>,
    max_n: Option<f64>,
) -> bool {
    if let Some(min) = min_len {
        if r.len() < min { return true; }
    }
    if let Some(max) = max_len {
        if r.len() > max { return true; }
    }
    if let Some(n) = max_n {
        let n_count = r.sequence.chars().filter(|&c| c == 'N' || c == 'n').count();
        let threshold = if n < 1.0 { (n * r.len() as f64) as usize } else { n as usize };
        if n_count > threshold { return true; }
    }
    false
}

fn write_record(
    writer: &mut dyn std::io::Write,
    record: &cutadapt_core::record::SequenceRecord,
    format: FileFormat,
) -> std::io::Result<()> {
    match format {
        FileFormat::Fastq => {
            let quals = record.qualities.as_deref().unwrap_or("");
            write!(writer, "@{}\n{}\n+\n{}\n", record.name, record.sequence, quals)
        }
        FileFormat::Fasta => {
            write!(writer, ">{}\n{}\n", record.name, record.sequence)
        }
    }
}

fn parse_paired_length(s: Option<&str>) -> Result<(Option<usize>, Option<usize>)> {
    match s {
        None => Ok((None, None)),
        Some(s) => {
            let vals = parse_lengths(s)?;
            let v1 = vals.first().and_then(|v| *v).map(|v| v as usize);
            let v2 = vals.get(1).and_then(|v| *v).map(|v| v as usize);
            Ok((v1, v2.or(v1)))
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn print_report(
    paired: bool,
    n_reads: u64,
    n_written: u64,
    bp_in1: u64,
    bp_in2: u64,
    bp_out1: u64,
    bp_out2: u64,
    reads_with_adapter1: u64,
    reads_with_adapter2: u64,
) {
    println!("=== Summary ===");
    println!();
    if paired {
        println!("Total read pairs processed:    {:>12}", fmt_num(n_reads));
        println!("  Read 1 with adapter:         {:>12} ({:.1}%)",
            fmt_num(reads_with_adapter1),
            pct(reads_with_adapter1, n_reads));
        println!("  Read 2 with adapter:         {:>12} ({:.1}%)",
            fmt_num(reads_with_adapter2),
            pct(reads_with_adapter2, n_reads));
        println!("Pairs written (passing filters): {:>10} ({:.1}%)",
            fmt_num(n_written), pct(n_written, n_reads));
        println!();
        let bp_in = bp_in1 + bp_in2;
        let bp_out = bp_out1 + bp_out2;
        println!("Total basepairs processed: {:>14} bp", fmt_num(bp_in));
        println!("  Read 1: {:>24} bp", fmt_num(bp_in1));
        println!("  Read 2: {:>24} bp", fmt_num(bp_in2));
        println!("Total written (filtered):  {:>14} bp ({:.1}%)",
            fmt_num(bp_out), pct(bp_out, bp_in));
        println!("  Read 1: {:>24} bp", fmt_num(bp_out1));
        println!("  Read 2: {:>24} bp", fmt_num(bp_out2));
    } else {
        println!("Total reads processed:         {:>12}", fmt_num(n_reads));
        println!("Reads with adapters:           {:>12} ({:.1}%)",
            fmt_num(reads_with_adapter1),
            pct(reads_with_adapter1, n_reads));
        println!("Reads written (passing filters): {:>10} ({:.1}%)",
            fmt_num(n_written), pct(n_written, n_reads));
        println!();
        println!("Total basepairs processed: {:>14} bp", fmt_num(bp_in1));
        println!("Total written (filtered):  {:>14} bp ({:.1}%)",
            fmt_num(bp_out1), pct(bp_out1, bp_in1));
    }
}

fn fmt_num(n: u64) -> String {
    let s = n.to_string();
    let mut result = String::new();
    for (i, c) in s.chars().rev().enumerate() {
        if i > 0 && i % 3 == 0 { result.push(','); }
        result.push(c);
    }
    result.chars().rev().collect()
}

fn pct(num: u64, den: u64) -> f64 {
    if den == 0 { 0.0 } else { num as f64 / den as f64 * 100.0 }
}

>>>>>>> 0ad2e08 (Wire up read processing pipeline in main.rs)
// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_cutoffs_single() {
        assert_eq!(parse_cutoffs("5").unwrap(), (0, 5));
    }

    #[test]
    fn test_parse_cutoffs_pair() {
        assert_eq!(parse_cutoffs("6,7").unwrap(), (6, 7));
    }

    #[test]
    fn test_parse_cutoffs_too_many() {
        assert!(parse_cutoffs("1,2,3").is_err());
    }

    #[test]
    fn test_parse_cutoffs_invalid() {
        assert!(parse_cutoffs("abc").is_err());
    }

    #[test]
    fn test_parse_lengths_single() {
        assert_eq!(parse_lengths("25").unwrap(), vec![Some(25)]);
    }

    #[test]
    fn test_parse_lengths_pair() {
        assert_eq!(parse_lengths("17:25").unwrap(), vec![Some(17), Some(25)]);
    }

    #[test]
    fn test_parse_lengths_trailing_colon() {
        assert_eq!(parse_lengths("25:").unwrap(), vec![Some(25), None]);
    }

    #[test]
    fn test_parse_lengths_leading_colon() {
        assert_eq!(parse_lengths(":25").unwrap(), vec![None, Some(25)]);
    }

    #[test]
    fn test_parse_lengths_double_colon_empty() {
        assert!(parse_lengths(":").is_err());
    }

    #[test]
    fn test_parse_lengths_too_many_colons() {
        assert!(parse_lengths("1:2:3").is_err());
    }
}
