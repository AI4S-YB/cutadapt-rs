//! File I/O for FASTA/FASTQ sequence files.
//!
//! Provides format detection, reading, and writing of sequence records,
//! with a structure that supports future compression (gzip, xz, bzip2).
//!
//! Ported from cutadapt/files.py.

use std::fmt;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Read, Write};
use std::collections::BTreeMap;
use std::sync::{mpsc, Arc, Mutex};
use std::thread::{self, JoinHandle};

use bzip2::read::BzDecoder;
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use xz2::read::XzDecoder;

use crate::record::SequenceRecord;

const DEFAULT_BUFFER_SIZE: usize = 1 << 20;

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
/// For compressed files (.gz, .bz2, .xz), decompresses the header first.
pub fn detect_file_format_from_path(path: &str) -> Result<FileFormat, FileFormatError> {
    if path == "-" {
        // Cannot peek at stdin without consuming bytes; default to FASTQ.
        return Ok(FileFormat::Fastq);
    }
    let opener = FileOpener::default();
    let mut reader = opener.open_read(path).map_err(FileFormatError::Io)?;
    let mut header = [0u8; 4];
    let n = reader.read(&mut header).map_err(FileFormatError::Io)?;
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

        records.push(SequenceRecord::from_parts(name, sequence, Some(qualities)));
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
                records.push(SequenceRecord::from_parts(
                    name,
                    std::mem::take(&mut current_seq),
                    None,
                ));
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
        records.push(SequenceRecord::from_parts(name, current_seq, None));
    }

    Ok(records)
}

fn trim_line_ending(line: &mut String) {
    if line.ends_with('\n') {
        line.pop();
        if line.ends_with('\r') {
            line.pop();
        }
    }
}

/// Parse up to `max_records` FASTQ records from a buffered reader.
fn read_fastq_batch(
    reader: &mut dyn BufRead,
    max_records: usize,
) -> Result<Vec<SequenceRecord>, io::Error> {
    let mut records = Vec::with_capacity(max_records);
    let mut header = String::new();
    let mut sequence = String::new();
    let mut separator = String::new();
    let mut qualities = String::new();

    while records.len() < max_records {
        header.clear();
        let n = reader.read_line(&mut header)?;
        if n == 0 {
            break;
        }
        trim_line_ending(&mut header);
        if header.is_empty() {
            continue;
        }
        if !header.starts_with('@') {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("Expected FASTQ header starting with '@', got: {header}"),
            ));
        }

        sequence.clear();
        if reader.read_line(&mut sequence)? == 0 {
            return Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                "Truncated FASTQ: missing sequence",
            ));
        }
        trim_line_ending(&mut sequence);

        separator.clear();
        if reader.read_line(&mut separator)? == 0 {
            return Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                "Truncated FASTQ: missing '+' line",
            ));
        }
        trim_line_ending(&mut separator);
        if !separator.starts_with('+') {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("Expected '+' separator in FASTQ, got: {separator}"),
            ));
        }

        qualities.clear();
        if reader.read_line(&mut qualities)? == 0 {
            return Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                "Truncated FASTQ: missing quality line",
            ));
        }
        trim_line_ending(&mut qualities);

        let seq_capacity = sequence.capacity();
        let qual_capacity = qualities.capacity();
        let sequence_owned = std::mem::take(&mut sequence);
        let qualities_owned = std::mem::take(&mut qualities);
        sequence = String::with_capacity(seq_capacity.max(sequence_owned.len()));
        qualities = String::with_capacity(qual_capacity.max(qualities_owned.len()));

        records.push(SequenceRecord::from_parts(
            header[1..].to_string(),
            sequence_owned,
            Some(qualities_owned),
        ));
    }

    Ok(records)
}

#[derive(Default)]
struct FastaReadState {
    current_name: Option<String>,
    current_sequence: String,
    reached_eof: bool,
}

/// Parse up to `max_records` FASTA records from a buffered reader, preserving
/// parser state across calls so the caller can stream batches.
fn read_fasta_batch(
    reader: &mut dyn BufRead,
    state: &mut FastaReadState,
    max_records: usize,
) -> Result<Vec<SequenceRecord>, io::Error> {
    let mut records = Vec::with_capacity(max_records);
    let mut line = String::new();

    while records.len() < max_records {
        line.clear();
        let n = reader.read_line(&mut line)?;
        if n == 0 {
            state.reached_eof = true;
            break;
        }
        trim_line_ending(&mut line);
        if line.is_empty() {
            continue;
        }
        if let Some(rest) = line.strip_prefix('>') {
            if let Some(name) = state.current_name.replace(rest.to_string()) {
                let sequence = std::mem::take(&mut state.current_sequence);
                records.push(SequenceRecord::from_parts(name, sequence, None));
                if records.len() == max_records {
                    return Ok(records);
                }
            }
        } else if line.starts_with('#') {
            continue;
        } else {
            state.current_sequence.push_str(&line);
        }
    }

    if state.reached_eof && records.len() < max_records {
        if let Some(name) = state.current_name.take() {
            let sequence = std::mem::take(&mut state.current_sequence);
            records.push(SequenceRecord::from_parts(name, sequence, None));
        }
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
    /// Buffer size used by buffered readers/writers.
    pub buffer_size: usize,
}

impl Default for FileOpener {
    fn default() -> Self {
        Self {
            compression_level: 1,
            threads: None,
            buffer_size: DEFAULT_BUFFER_SIZE,
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

struct GzipJob {
    sequence: usize,
    data: Vec<u8>,
}

struct ParallelGzipWriter {
    sink: BufWriter<File>,
    sender: Option<mpsc::Sender<Option<GzipJob>>>,
    receiver: mpsc::Receiver<(usize, io::Result<Vec<u8>>)>,
    workers: Vec<JoinHandle<()>>,
    pending: BTreeMap<usize, Vec<u8>>,
    buffer: Vec<u8>,
    chunk_size: usize,
    max_in_flight: usize,
    next_sequence: usize,
    next_to_write: usize,
    in_flight: usize,
}

impl ParallelGzipWriter {
    fn new(
        file: File,
        compression_level: u32,
        worker_count: usize,
        buffer_size: usize,
    ) -> io::Result<Self> {
        let chunk_size = buffer_size.max(DEFAULT_BUFFER_SIZE * 8);
        let sink = BufWriter::with_capacity(buffer_size, file);
        let (sender, job_receiver) = mpsc::channel::<Option<GzipJob>>();
        let (result_sender, receiver) = mpsc::channel::<(usize, io::Result<Vec<u8>>)>();
        let shared_receiver = Arc::new(Mutex::new(job_receiver));
        let mut workers = Vec::with_capacity(worker_count);

        for _ in 0..worker_count {
            let rx = Arc::clone(&shared_receiver);
            let tx = result_sender.clone();
            workers.push(thread::spawn(move || {
                loop {
                    let job = match rx.lock().unwrap().recv() {
                        Ok(job) => job,
                        Err(_) => break,
                    };
                    let Some(job) = job else {
                        break;
                    };
                    let result = compress_gzip_member(job.data, compression_level);
                    if tx.send((job.sequence, result)).is_err() {
                        break;
                    }
                }
            }));
        }
        drop(result_sender);

        Ok(Self {
            sink,
            sender: Some(sender),
            receiver,
            workers,
            pending: BTreeMap::new(),
            buffer: Vec::with_capacity(chunk_size),
            chunk_size,
            max_in_flight: worker_count.max(1),
            next_sequence: 0,
            next_to_write: 0,
            in_flight: 0,
        })
    }

    fn enqueue_chunk(&mut self, data: Vec<u8>) -> io::Result<()> {
        if data.is_empty() {
            return Ok(());
        }

        let sender = self.sender.as_ref().ok_or_else(|| {
            io::Error::new(io::ErrorKind::BrokenPipe, "gzip writer has been closed")
        })?;
        sender
            .send(Some(GzipJob {
                sequence: self.next_sequence,
                data,
            }))
            .map_err(|_| io::Error::new(io::ErrorKind::BrokenPipe, "gzip worker channel closed"))?;
        self.next_sequence += 1;
        self.in_flight += 1;

        if self.in_flight >= self.max_in_flight {
            self.collect_one_blocking()?;
        }
        self.collect_available()?;
        Ok(())
    }

    fn flush_pending_buffer(&mut self) -> io::Result<()> {
        if self.buffer.is_empty() {
            return Ok(());
        }
        let payload = std::mem::take(&mut self.buffer);
        self.enqueue_chunk(payload)
    }

    fn process_result(&mut self, sequence: usize, result: io::Result<Vec<u8>>) -> io::Result<()> {
        self.in_flight = self.in_flight.saturating_sub(1);
        let data = result?;
        self.pending.insert(sequence, data);

        while let Some(chunk) = self.pending.remove(&self.next_to_write) {
            self.sink.write_all(&chunk)?;
            self.next_to_write += 1;
        }
        Ok(())
    }

    fn collect_available(&mut self) -> io::Result<()> {
        while let Ok((sequence, result)) = self.receiver.try_recv() {
            self.process_result(sequence, result)?;
        }
        Ok(())
    }

    fn collect_one_blocking(&mut self) -> io::Result<()> {
        let (sequence, result) = self
            .receiver
            .recv()
            .map_err(|_| io::Error::new(io::ErrorKind::BrokenPipe, "gzip worker channel closed"))?;
        self.process_result(sequence, result)
    }

    fn finish_open_chunks(&mut self) -> io::Result<()> {
        self.flush_pending_buffer()?;
        while self.in_flight > 0 {
            self.collect_one_blocking()?;
        }
        self.collect_available()?;
        self.sink.flush()
    }

    fn shutdown_workers(&mut self) {
        if let Some(sender) = self.sender.take() {
            for _ in 0..self.workers.len() {
                let _ = sender.send(None);
            }
        }
        for worker in self.workers.drain(..) {
            let _ = worker.join();
        }
    }
}

impl Write for ParallelGzipWriter {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        if buf.is_empty() {
            return Ok(0);
        }

        if self.buffer.is_empty() && buf.len() >= self.chunk_size {
            self.enqueue_chunk(buf.to_vec())?;
            return Ok(buf.len());
        }

        self.buffer.extend_from_slice(buf);
        while self.buffer.len() >= self.chunk_size {
            let tail = self.buffer.split_off(self.chunk_size);
            let payload = std::mem::replace(&mut self.buffer, tail);
            self.enqueue_chunk(payload)?;
        }
        Ok(buf.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        self.finish_open_chunks()
    }
}

impl Drop for ParallelGzipWriter {
    fn drop(&mut self) {
        let _ = self.finish_open_chunks();
        self.shutdown_workers();
    }
}

fn compress_gzip_member(data: Vec<u8>, compression_level: u32) -> io::Result<Vec<u8>> {
    let mut encoder = GzEncoder::new(Vec::new(), flate2::Compression::new(compression_level));
    encoder.write_all(&data)?;
    encoder.finish()
}

impl FileOpener {
    pub fn new(compression_level: i32, threads: Option<usize>) -> Self {
        Self {
            compression_level,
            threads,
            buffer_size: DEFAULT_BUFFER_SIZE,
        }
    }

    pub fn with_buffer_size(mut self, buffer_size: usize) -> Self {
        self.buffer_size = buffer_size.max(8 * 1024);
        self
    }

    fn effective_buffer_size(&self) -> usize {
        self.buffer_size.max(8 * 1024)
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
    /// Supports `"-"` for stdin. Compressed formats (.gz, .bz2, .xz) are
    /// auto-detected from the file extension.
    pub fn open_read(&self, path: &str) -> io::Result<Box<dyn BufRead>> {
        let capacity = self.effective_buffer_size();
        if path == "-" {
            return Ok(Box::new(BufReader::with_capacity(capacity, io::stdin())));
        }
        let file = File::open(path)?;
        match Self::detect_compression(path) {
            Compression::None => Ok(Box::new(BufReader::with_capacity(capacity, file))),
            Compression::Gzip => Ok(Box::new(BufReader::with_capacity(
                capacity,
                MultiGzDecoder::new(BufReader::with_capacity(capacity, file)),
            ))),
            Compression::Bzip2 => Ok(Box::new(BufReader::with_capacity(
                capacity,
                BzDecoder::new(BufReader::with_capacity(capacity, file)),
            ))),
            Compression::Xz => Ok(Box::new(BufReader::with_capacity(
                capacity,
                XzDecoder::new(BufReader::with_capacity(capacity, file)),
            ))),
        }
    }

    /// Open a file for writing. Returns a boxed `Write`.
    ///
    /// Supports `"-"` for stdout. Gzip compression is applied for `.gz` files.
    pub fn open_write(&self, path: &str) -> io::Result<Box<dyn Write>> {
        let capacity = self.effective_buffer_size();
        if path == "-" {
            return Ok(Box::new(BufWriter::with_capacity(capacity, io::stdout())));
        }
        match Self::detect_compression(path) {
            Compression::None => {
                let file = File::create(path)?;
                Ok(Box::new(BufWriter::with_capacity(capacity, file)))
            }
            Compression::Gzip => {
                let level = flate2::Compression::new(self.compression_level.max(1) as u32);
                let worker_count = self.threads.unwrap_or(1);
                if worker_count > 1 {
                    Ok(Box::new(ParallelGzipWriter::new(
                        File::create(path)?,
                        self.compression_level.max(1) as u32,
                        worker_count,
                        capacity,
                    )?))
                } else {
                    let file = BufWriter::with_capacity(capacity, File::create(path)?);
                    let encoder = GzEncoder::new(file, level);
                    Ok(Box::new(BufWriter::with_capacity(capacity, encoder)))
                }
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
        self.open_with_opener(&FileOpener::default())
    }

    pub fn open_with_opener(&self, opener: &FileOpener) -> io::Result<InputFiles> {
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
            fasta_states: std::iter::repeat_with(FastaReadState::default)
                .take(self.paths.len())
                .collect(),
        })
    }
}

/// Opened input files ready for reading.
pub struct InputFiles {
    pub readers: Vec<Box<dyn BufRead>>,
    pub format: FileFormat,
    pub interleaved: bool,
    fasta_states: Vec<FastaReadState>,
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

    pub fn read_batch(
        &mut self,
        reader_index: usize,
        batch_size: usize,
    ) -> io::Result<Vec<SequenceRecord>> {
        let reader = self
            .readers
            .get_mut(reader_index)
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "Reader index out of range"))?;
        match self.format {
            FileFormat::Fastq => read_fastq_batch(reader.as_mut(), batch_size),
            FileFormat::Fasta => {
                let state = self
                    .fasta_states
                    .get_mut(reader_index)
                    .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "Reader index out of range"))?;
                read_fasta_batch(reader.as_mut(), state, batch_size)
            }
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

    #[test]
    fn test_parallel_gzip_writer_roundtrip() {
        let dir = std::env::temp_dir();
        let path = dir.join("cutadapt_test_parallel_output.fastq.gz");
        let path_str = path.to_str().unwrap();

        let opener = FileOpener::new(1, Some(2)).with_buffer_size(1024);
        {
            let mut writer = opener.open_write(path_str).unwrap();
            let rec1 = SequenceRecord::new("test1", "ACGT", Some("IIII"));
            let rec2 = SequenceRecord::new("test2", "TGCA", Some("HHHH"));
            write_fastq(&mut writer, &rec1).unwrap();
            write_fastq(&mut writer, &rec2).unwrap();
            writer.flush().unwrap();
        }

        let mut reader = opener.open_read(path_str).unwrap();
        let parsed = read_fastq(&mut reader).unwrap();
        assert_eq!(parsed.len(), 2);
        assert_eq!(parsed[0].name, "test1");
        assert_eq!(parsed[0].sequence, "ACGT");
        assert_eq!(parsed[1].name, "test2");
        assert_eq!(parsed[1].qualities.as_deref(), Some("HHHH"));

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
