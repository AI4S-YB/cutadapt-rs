//! File I/O for FASTA/FASTQ sequence files.
//!
//! Provides format detection, reading, and writing of sequence records,
//! with a structure that supports future compression (gzip, xz, bzip2).
//!
//! Ported from cutadapt/files.py.

use std::fmt;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Read, Write};

use crate::record::SequenceRecord;

// ---------------------------------------------------------------------------
// FileFormat
// ---------------------------------------------------------------------------

/// Sequence file format.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FileFormat {
    Fasta,
    Fastq,
}

impl FileFormat {
    /// Returns `true` if this format carries per-base quality scores.
    pub fn has_qualities(&self) -> bool {
        matches!(self, FileFormat::Fastq)
    }
}

impl fmt::Display for FileFormat {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            FileFormat::Fasta => write!(f, "FASTA"),
            FileFormat::Fastq => write!(f, "FASTQ"),
        }
    }
}

// ---------------------------------------------------------------------------
// Format detection
// ---------------------------------------------------------------------------

/// Detect the file format from the first bytes of an already-opened reader.
///
/// Rules (matching the Python implementation):
/// - Starts with `@` or is empty → FASTQ
/// - Starts with `>` or `#` → FASTA
/// - Otherwise → error
pub fn detect_file_format(buf: &[u8]) -> Result<FileFormat, FileFormatError> {
    if buf.is_empty() || buf[0] == b'@' {
        // Pretend FASTQ for empty input (matches Python behaviour)
        Ok(FileFormat::Fastq)
    } else if buf[0] == b'>' || buf[0] == b'#' {
        Ok(FileFormat::Fasta)
    } else {
        Err(FileFormatError::Unknown {
            first_bytes: buf[..buf.len().min(4)].to_vec(),
        })
    }
}

/// Detect format by opening a file, reading a small header, then closing it.
pub fn detect_file_format_from_path(path: &str) -> Result<FileFormat, FileFormatError> {
    if path == "-" {
        // Cannot peek at stdin without consuming bytes; default to FASTQ.
        return Ok(FileFormat::Fastq);
    }
    let mut f = File::open(path).map_err(FileFormatError::Io)?;
    let mut header = [0u8; 4];
    let n = f.read(&mut header).map_err(FileFormatError::Io)?;
    detect_file_format(&header[..n])
}

/// Errors that can occur during format detection.
#[derive(Debug)]
pub enum FileFormatError {
    Unknown { first_bytes: Vec<u8> },
    Io(io::Error),
}

impl fmt::Display for FileFormatError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            FileFormatError::Unknown { first_bytes } => {
                write!(
                    f,
                    "Input file format not recognized. The file starts with {:?}, \
                     but files in supported formats start with '>' (FASTA) or '@' (FASTQ)",
                    first_bytes
                )
            }
            FileFormatError::Io(e) => write!(f, "I/O error during format detection: {e}"),
        }
    }
}

impl std::error::Error for FileFormatError {}

// ---------------------------------------------------------------------------
// Reading helpers
// ---------------------------------------------------------------------------

/// Parse all FASTQ records from a buffered reader.
///
/// FASTQ format (4 lines per record):
/// ```text
/// @name
/// SEQUENCE
/// +
/// QUALITIES
/// ```
pub fn read_fastq(reader: impl BufRead) -> Result<Vec<SequenceRecord>, io::Error> {
    let mut records = Vec::new();
    let mut lines = reader.lines();

    loop {
        // Line 1: header starting with '@'
        let header = match lines.next() {
            Some(Ok(line)) => line,
            Some(Err(e)) => return Err(e),
            None => break, // EOF
        };
        if header.is_empty() {
            // Skip blank trailing lines
            continue;
        }
        if !header.starts_with('@') {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("Expected FASTQ header starting with '@', got: {header}"),
            ));
        }
        let name = header[1..].to_string();

        // Line 2: sequence
        let sequence = lines
            .next()
            .ok_or_else(|| {
                io::Error::new(io::ErrorKind::UnexpectedEof, "Truncated FASTQ: missing sequence")
            })?
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

        // Line 3: '+' separator
        let sep = lines
            .next()
            .ok_or_else(|| {
                io::Error::new(io::ErrorKind::UnexpectedEof, "Truncated FASTQ: missing '+' line")
            })?
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        if !sep.starts_with('+') {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("Expected '+' separator in FASTQ, got: {sep}"),
            ));
        }

        // Line 4: quality scores
        let qualities = lines
            .next()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::UnexpectedEof,
                    "Truncated FASTQ: missing quality line",
                )
            })?
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

        records.push(SequenceRecord::new(&name, &sequence, Some(&qualities)));
    }

    Ok(records)
}

/// Parse all FASTA records from a buffered reader.
///
/// FASTA format:
/// ```text
/// >name optional description
/// SEQUENCE
/// CONTINUED...
/// ```
/// Multi-line sequences are concatenated.
pub fn read_fasta(reader: impl BufRead) -> Result<Vec<SequenceRecord>, io::Error> {
    let mut records = Vec::new();
    let mut current_name: Option<String> = None;
    let mut current_seq = String::new();

    for line_result in reader.lines() {
        let line = line_result?;
        if line.is_empty() {
            continue;
        }
        if line.starts_with('>') {
            // Flush previous record
            if let Some(name) = current_name.take() {
                records.push(SequenceRecord::new(&name, &current_seq, None));
                current_seq.clear();
            }
            current_name = Some(line[1..].to_string());
        } else if line.starts_with('#') {
            // Comment line — skip
            continue;
        } else {
            current_seq.push_str(&line);
        }
    }
    // Flush last record
    if let Some(name) = current_name.take() {
        records.push(SequenceRecord::new(&name, &current_seq, None));
    }

    Ok(records)
}

// ---------------------------------------------------------------------------
// Writing helpers
// ---------------------------------------------------------------------------

/// Write a single record in FASTQ format.
pub fn write_fastq(writer: &mut impl Write, record: &SequenceRecord) -> io::Result<()> {
    let quals = record.qualities.as_deref().unwrap_or("");
    write!(writer, "@{}\n{}\n+\n{}\n", record.name, record.sequence, quals)
}

/// Write a single record in FASTA format.
pub fn write_fasta(writer: &mut impl Write, record: &SequenceRecord) -> io::Result<()> {
    write!(writer, ">{}\n{}\n", record.name, record.sequence)
}

// ---------------------------------------------------------------------------
// FileOpener
// ---------------------------------------------------------------------------

/// Configures how files are opened (compression level, thread count).
///
/// Currently only handles plain (uncompressed) files.  Compressed formats
/// (.gz, .xz, .bz2) are detected by extension so the structure is ready
/// for future wiring of flate2 / xz2 / bzip2 crates.
#[derive(Debug, Clone)]
pub struct FileOpener {
    /// Compression level for writing (1-9). Ignored for uncompressed files.
    pub compression_level: i32,
    /// Number of threads for external compression.
    /// `None` means use a sensible default.
    pub threads: Option<usize>,
}

impl Default for FileOpener {
    fn default() -> Self {
        Self {
            compression_level: 1,
            threads: None,
        }
    }
}

/// Compression scheme inferred from a file extension.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Compression {
    None,
    Gzip,
    Bzip2,
    Xz,
}

impl FileOpener {
    pub fn new(compression_level: i32, threads: Option<usize>) -> Self {
        Self {
            compression_level,
            threads,
        }
    }

    /// Infer compression from the file extension.
    pub fn detect_compression(path: &str) -> Compression {
        if path.ends_with(".gz") {
            Compression::Gzip
        } else if path.ends_with(".bz2") {
            Compression::Bzip2
        } else if path.ends_with(".xz") {
            Compression::Xz
        } else {
            Compression::None
        }
    }

    /// Open a file for reading. Returns a boxed `BufRead`.
    ///
    /// Supports `"-"` for stdin. Compressed formats are detected but
    /// currently only plain files are handled (compressed support is
    /// stubbed and will error).
    pub fn open_read(&self, path: &str) -> io::Result<Box<dyn BufRead>> {
        if path == "-" {
            return Ok(Box::new(BufReader::new(io::stdin())));
        }
        let compression = Self::detect_compression(path);
        match compression {
            Compression::None => {
                let file = File::open(path)?;
                Ok(Box::new(BufReader::new(file)))
            }
            other => Err(io::Error::new(
                io::ErrorKind::Unsupported,
                format!("Compression {:?} not yet implemented for reading", other),
            )),
        }
    }

    /// Open a file for writing. Returns a boxed `Write`.
    ///
    /// Supports `"-"` for stdout. Compressed formats are detected but
    /// currently only plain files are handled.
    pub fn open_write(&self, path: &str) -> io::Result<Box<dyn Write>> {
        if path == "-" {
            return Ok(Box::new(BufWriter::new(io::stdout())));
        }
        let compression = Self::detect_compression(path);
        match compression {
            Compression::None => {
                let file = File::create(path)?;
                Ok(Box::new(BufWriter::new(file)))
            }
            other => Err(io::Error::new(
                io::ErrorKind::Unsupported,
                format!("Compression {:?} not yet implemented for writing", other),
            )),
        }
    }
}

// ---------------------------------------------------------------------------
// InputPaths / InputFiles
// ---------------------------------------------------------------------------

/// Holds the paths to input files (before they are opened).
#[derive(Debug, Clone)]
pub struct InputPaths {
    pub paths: Vec<String>,
    pub interleaved: bool,
}

impl InputPaths {
    pub fn new(paths: Vec<String>, interleaved: bool) -> Self {
        Self { paths, interleaved }
    }

    /// Open all input paths and return an `InputFiles`.
    pub fn open(&self) -> io::Result<InputFiles> {
        let opener = FileOpener::default();
        let mut readers: Vec<Box<dyn BufRead>> = Vec::with_capacity(self.paths.len());
        let mut format: Option<FileFormat> = None;

        for path in &self.paths {
            let reader = opener.open_read(path)?;
            // Detect format from the first file if we haven't yet
            if format.is_none() && path != "-" {
                format = detect_file_format_from_path(path).ok();
            }
            readers.push(reader);
        }

        Ok(InputFiles {
            readers,
            format: format.unwrap_or(FileFormat::Fastq),
            interleaved: self.interleaved,
        })
    }
}

/// Opened input files ready for reading.
pub struct InputFiles {
    pub readers: Vec<Box<dyn BufRead>>,
    pub format: FileFormat,
    pub interleaved: bool,
}

impl InputFiles {
    /// Read all records from the first reader.
    ///
    /// For paired-end data the caller should also process additional readers;
    /// this convenience method handles the common single-file case.
    pub fn read_all(&mut self) -> io::Result<Vec<SequenceRecord>> {
        if self.readers.is_empty() {
            return Ok(Vec::new());
        }
        let reader = &mut self.readers[0];
        match self.format {
            FileFormat::Fastq => read_fastq(reader),
            FileFormat::Fasta => read_fasta(reader),
        }
    }
}

// ---------------------------------------------------------------------------
// OutputFiles
// ---------------------------------------------------------------------------

/// Manages output files for writing sequence records.
pub struct OutputFiles {
    writers: Vec<Box<dyn Write>>,
    format: FileFormat,
    interleaved: bool,
    opener: FileOpener,
}

impl OutputFiles {
    pub fn new(qualities: bool, interleaved: bool, file_opener: Option<FileOpener>) -> Self {
        let format = if qualities {
            FileFormat::Fastq
        } else {
            FileFormat::Fasta
        };
        Self {
            writers: Vec::new(),
            format,
            interleaved,
            opener: file_opener.unwrap_or_default(),
        }
    }

    /// Open a file path for record writing and add it to the internal list.
    /// Returns the index of the newly opened writer.
    pub fn open_record_writer(&mut self, path: &str) -> io::Result<usize> {
        let writer = self.opener.open_write(path)?;
        let idx = self.writers.len();
        self.writers.push(writer);
        Ok(idx)
    }

    /// Write a record to the writer at the given index.
    pub fn write_record(&mut self, index: usize, record: &SequenceRecord) -> io::Result<()> {
        let writer = self
            .writers
            .get_mut(index)
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "Writer index out of range"))?;
        match self.format {
            FileFormat::Fastq => write_fastq(writer, record),
            FileFormat::Fasta => write_fasta(writer, record),
        }
    }

    /// Write a record using format auto-detection: FASTQ if qualities are
    /// present, FASTA otherwise. Uses the writer at the given index.
    pub fn write_record_auto(&mut self, index: usize, record: &SequenceRecord) -> io::Result<()> {
        let writer = self
            .writers
            .get_mut(index)
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "Writer index out of range"))?;
        if record.qualities.is_some() {
            write_fastq(writer, record)
        } else {
            write_fasta(writer, record)
        }
    }

    /// Flush and return the format being used.
    pub fn format(&self) -> FileFormat {
        self.format
    }

    /// Whether output is interleaved.
    pub fn interleaved(&self) -> bool {
        self.interleaved
    }

    /// Flush all writers.
    pub fn flush(&mut self) -> io::Result<()> {
        for w in &mut self.writers {
            w.flush()?;
        }
        Ok(())
    }

    /// Close (flush and drop) all writers.
    pub fn close(mut self) -> io::Result<()> {
        self.flush()?;
        // Dropping `self` closes file handles.
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    // -- FileFormat ---------------------------------------------------------

    #[test]
    fn test_file_format_has_qualities() {
        assert!(FileFormat::Fastq.has_qualities());
        assert!(!FileFormat::Fasta.has_qualities());
    }

    #[test]
    fn test_file_format_display() {
        assert_eq!(format!("{}", FileFormat::Fasta), "FASTA");
        assert_eq!(format!("{}", FileFormat::Fastq), "FASTQ");
    }

    // -- detect_file_format -------------------------------------------------

    #[test]
    fn test_detect_fastq() {
        assert_eq!(detect_file_format(b"@SEQ_ID").unwrap(), FileFormat::Fastq);
    }

    #[test]
    fn test_detect_fasta_gt() {
        assert_eq!(detect_file_format(b">SEQ_ID").unwrap(), FileFormat::Fasta);
    }

    #[test]
    fn test_detect_fasta_comment() {
        assert_eq!(detect_file_format(b"# comment").unwrap(), FileFormat::Fasta);
    }

    #[test]
    fn test_detect_empty_defaults_to_fastq() {
        assert_eq!(detect_file_format(b"").unwrap(), FileFormat::Fastq);
    }

    #[test]
    fn test_detect_unknown() {
        assert!(detect_file_format(b"XYZ").is_err());
    }

    // -- read_fastq ---------------------------------------------------------

    #[test]
    fn test_read_fastq_single() {
        let data = b"@read1\nACGT\n+\nIIII\n";
        let records = read_fastq(Cursor::new(data)).unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].name, "read1");
        assert_eq!(records[0].sequence, "ACGT");
        assert_eq!(records[0].qualities.as_deref(), Some("IIII"));
    }

    #[test]
    fn test_read_fastq_multiple() {
        let data = b"@r1\nACGT\n+\nIIII\n@r2\nTGCA\n+\nHHHH\n";
        let records = read_fastq(Cursor::new(data)).unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].name, "r1");
        assert_eq!(records[1].name, "r2");
        assert_eq!(records[1].sequence, "TGCA");
    }

    #[test]
    fn test_read_fastq_empty() {
        let data = b"";
        let records = read_fastq(Cursor::new(data)).unwrap();
        assert!(records.is_empty());
    }

    #[test]
    fn test_read_fastq_bad_header() {
        let data = b"read1\nACGT\n+\nIIII\n";
        let result = read_fastq(Cursor::new(data));
        assert!(result.is_err());
    }

    #[test]
    fn test_read_fastq_truncated() {
        let data = b"@read1\nACGT\n";
        let result = read_fastq(Cursor::new(data));
        assert!(result.is_err());
    }

    // -- read_fasta ---------------------------------------------------------

    #[test]
    fn test_read_fasta_single() {
        let data = b">read1\nACGT\n";
        let records = read_fasta(Cursor::new(data)).unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].name, "read1");
        assert_eq!(records[0].sequence, "ACGT");
        assert!(records[0].qualities.is_none());
    }

    #[test]
    fn test_read_fasta_multiline() {
        let data = b">read1\nACGT\nTGCA\n";
        let records = read_fasta(Cursor::new(data)).unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].sequence, "ACGTTGCA");
    }

    #[test]
    fn test_read_fasta_multiple() {
        let data = b">r1\nACGT\n>r2\nTGCA\n";
        let records = read_fasta(Cursor::new(data)).unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].name, "r1");
        assert_eq!(records[0].sequence, "ACGT");
        assert_eq!(records[1].name, "r2");
        assert_eq!(records[1].sequence, "TGCA");
    }

    #[test]
    fn test_read_fasta_with_comments() {
        let data = b"# comment\n>r1\nACGT\n";
        let records = read_fasta(Cursor::new(data)).unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].name, "r1");
    }

    #[test]
    fn test_read_fasta_empty() {
        let data = b"";
        let records = read_fasta(Cursor::new(data)).unwrap();
        assert!(records.is_empty());
    }

    // -- write_fastq / write_fasta ------------------------------------------

    #[test]
    fn test_write_fastq() {
        let rec = SequenceRecord::new("read1", "ACGT", Some("IIII"));
        let mut buf = Vec::new();
        write_fastq(&mut buf, &rec).unwrap();
        assert_eq!(buf, b"@read1\nACGT\n+\nIIII\n");
    }

    #[test]
    fn test_write_fasta() {
        let rec = SequenceRecord::new("read1", "ACGT", None);
        let mut buf = Vec::new();
        write_fasta(&mut buf, &rec).unwrap();
        assert_eq!(buf, b">read1\nACGT\n");
    }

    // -- round-trip ---------------------------------------------------------

    #[test]
    fn test_fastq_roundtrip() {
        let original = vec![
            SequenceRecord::new("r1", "ACGTACGT", Some("IIIIIIII")),
            SequenceRecord::new("r2 some description", "TGCA", Some("HHHH")),
        ];
        let mut buf = Vec::new();
        for rec in &original {
            write_fastq(&mut buf, rec).unwrap();
        }
        let parsed = read_fastq(Cursor::new(&buf)).unwrap();
        assert_eq!(parsed, original);
    }

    #[test]
    fn test_fasta_roundtrip() {
        let original = vec![
            SequenceRecord::new("r1", "ACGTACGT", None),
            SequenceRecord::new("r2 desc", "TGCA", None),
        ];
        let mut buf = Vec::new();
        for rec in &original {
            write_fasta(&mut buf, rec).unwrap();
        }
        let parsed = read_fasta(Cursor::new(&buf)).unwrap();
        assert_eq!(parsed, original);
    }

    // -- FileOpener ---------------------------------------------------------

    #[test]
    fn test_detect_compression() {
        assert_eq!(FileOpener::detect_compression("reads.fq"), Compression::None);
        assert_eq!(FileOpener::detect_compression("reads.fq.gz"), Compression::Gzip);
        assert_eq!(FileOpener::detect_compression("reads.fq.bz2"), Compression::Bzip2);
        assert_eq!(FileOpener::detect_compression("reads.fq.xz"), Compression::Xz);
    }

    #[test]
    fn test_file_opener_default() {
        let opener = FileOpener::default();
        assert_eq!(opener.compression_level, 1);
        assert_eq!(opener.threads, None);
    }

    // -- OutputFiles --------------------------------------------------------

    #[test]
    fn test_output_files_write_to_buffer() {
        // Use a temporary file to verify OutputFiles end-to-end.
        let dir = std::env::temp_dir();
        let path = dir.join("cutadapt_test_output.fq");
        let path_str = path.to_str().unwrap();

        let mut out = OutputFiles::new(true, false, None);
        let idx = out.open_record_writer(path_str).unwrap();
        let rec = SequenceRecord::new("testread", "ACGT", Some("IIII"));
        out.write_record(idx, &rec).unwrap();
        out.close().unwrap();

        // Read back
        let contents = std::fs::read_to_string(&path).unwrap();
        assert_eq!(contents, "@testread\nACGT\n+\nIIII\n");

        // Clean up
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn test_output_files_fasta() {
        let dir = std::env::temp_dir();
        let path = dir.join("cutadapt_test_output.fa");
        let path_str = path.to_str().unwrap();

        let mut out = OutputFiles::new(false, false, None);
        let idx = out.open_record_writer(path_str).unwrap();
        let rec = SequenceRecord::new("testread", "ACGT", None);
        out.write_record(idx, &rec).unwrap();
        out.close().unwrap();

        let contents = std::fs::read_to_string(&path).unwrap();
        assert_eq!(contents, ">testread\nACGT\n");

        let _ = std::fs::remove_file(&path);
    }

    // -- InputPaths with real file ------------------------------------------

    #[test]
    fn test_input_paths_fastq_file() {
        let dir = std::env::temp_dir();
        let path = dir.join("cutadapt_test_input.fq");
        std::fs::write(&path, "@r1\nACGT\n+\nIIII\n@r2\nTTTT\n+\nHHHH\n").unwrap();

        let input_paths = InputPaths::new(vec![path.to_str().unwrap().to_string()], false);
        let mut input_files = input_paths.open().unwrap();
        assert_eq!(input_files.format, FileFormat::Fastq);

        let records = input_files.read_all().unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].name, "r1");
        assert_eq!(records[1].sequence, "TTTT");

        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn test_input_paths_fasta_file() {
        let dir = std::env::temp_dir();
        let path = dir.join("cutadapt_test_input.fa");
        std::fs::write(&path, ">r1\nACGT\n>r2\nTGCA\n").unwrap();

        let input_paths = InputPaths::new(vec![path.to_str().unwrap().to_string()], false);
        let mut input_files = input_paths.open().unwrap();
        assert_eq!(input_files.format, FileFormat::Fasta);

        let records = input_files.read_all().unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].name, "r1");
        assert_eq!(records[1].sequence, "TGCA");

        let _ = std::fs::remove_file(&path);
    }
}
