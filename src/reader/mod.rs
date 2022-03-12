use std::io::Read;

#[macro_use]
mod errors;
mod nom_parsers;
mod streaming_parser;
use self::streaming_parser::StreamParser;
use crate::seq::{Location, Seq};

pub use crate::errors::GbParserError;

#[derive(Debug)]
pub struct SeqReader<T: Read> {
    parser: StreamParser<T>,
}

impl<T: Read> Iterator for SeqReader<T> {
    type Item = Result<Seq, GbParserError>;

    fn next(&mut self) -> Option<Result<Seq, GbParserError>> {
        match self.parser.read_one_record() {
            Ok(Some(seq)) => Some(Ok(seq)),
            Ok(None) => None,
            Err(e) => Some(Err(e)),
        }
    }
}

const READ_BUF_SIZE: usize = 64 * 1024;

impl<T: Read> SeqReader<T> {
    /// Parse a stream one `Seq` at a time
    pub fn new(data: T) -> SeqReader<T> {
        SeqReader {
            parser: StreamParser::new(data, READ_BUF_SIZE),
        }
    }
}

/// Convenience method to parse an entire file at once. Uses the streaming parser.
pub fn parse_file<P: AsRef<::std::path::Path>>(path: P) -> Result<Vec<Seq>, GbParserError> {
    let file = ::std::fs::File::open(path)?;
    SeqReader::new(file).collect()
}

/// Parse an entire genbank file provided as a slice. Might be slightly faster
/// than the streaming parser used by `parse_file` and `SeqReader::from_stream` 
/// since less copying of data is required, however not as well tested. I recommend using
/// `SeqReader` instead. I've mainly left this here for benchmarking purposes.
pub fn parse_slice(data: &[u8]) -> Result<Vec<Seq>, GbParserError> {
    let res = nom_parsers::gb_records(data);
    match res {
        Ok((_, o)) => Ok(o),
        Err(e) => {
            Err(GbParserError::SyntaxError(
                format!("{:?}", e),
            ))
        }
    }
}

/// used by `Location::from_gb_format`
pub (crate) fn parse_location(data: &[u8]) -> Result<Location, GbParserError> {
    let res = nom_parsers::location(nom::types::CompleteByteSlice(data));
    match res {
        Ok((_, o)) => Ok(o),
        Err(e) => {
            Err(GbParserError::SyntaxError(
                format!("{:?}", e),
            ))
        }
    }
}