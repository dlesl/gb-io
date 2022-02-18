use crate::reader::nom_parsers::{
    any_field, base_count, contig_text, double_slash, feature, features_header, fill_seq_fields,
    line_ending_type_hack, locus, origin_tag, skip_preamble,
};
use nom::{self, AsChar, IResult, Offset};
use std::cmp;
use std::io::Error as IoError;
use std::io::Read;
use std::io::Result as IoResult;

use crate::seq::*;

use crate::errors::GbParserError;

extern crate circular;

#[derive(Debug)]
pub struct StreamParser<T: Read> {
    buffer: circular::Buffer,
    stream: T,
    capacity: usize,
    is_eof: bool,
}

// We use this private error type rather than nom's errors, so that we can own
// the input slice to give "context" even once the slice we were parsing is gone

const MAX_CONTEXT_BYTES: usize = 50; // maximum length of the Vec in the StreamParser
                                     // variant to avoid cloning massive input
                                     // slices

enum StreamParserError {
    Io(IoError),
    StreamParser(Option<Vec<u8>>, nom::ErrorKind),
    EOF,
}

impl From<IoError> for StreamParserError {
    fn from(e: IoError) -> StreamParserError {
        StreamParserError::Io(e)
    }
}

impl From<StreamParserError> for GbParserError {
    fn from(e: StreamParserError) -> GbParserError {
        match e {
            StreamParserError::Io(e) => GbParserError::from(e),
            StreamParserError::EOF => GbParserError::SyntaxError("Unexpected EOF".into()),
            StreamParserError::StreamParser(Some(context), e) => {
                GbParserError::SyntaxError(format!(
                    "Error {:?} while parsing [{}]",
                    e,
                    String::from_utf8_lossy(&context)
                ))
            }
            StreamParserError::StreamParser(None, e) => {
                GbParserError::SyntaxError(format!("Parse error: {:?}", e))
            }
        }
    }
}

impl<T: Read> StreamParser<T> {
    pub fn new(stream: T, capacity: usize) -> StreamParser<T> {
        StreamParser {
            stream,
            capacity,
            buffer: circular::Buffer::with_capacity(capacity),
            is_eof: false,
        }
    }

    fn fill_buffer(&mut self) -> IoResult<usize> {
        if self.is_eof() {
            return Ok(0);
        }
        // if we're requesting a buffer refill when the buffer's full, we need
        // to grow it.
        if self.buffer.available_space() == 0 {
            self.capacity *= 2;
            self.buffer.grow(self.capacity);
            debug!("Increasing read buffer capacity to {} b", self.capacity);
        }
        let bytes_read = self.stream.read(self.buffer.space())?;
        if bytes_read == 0 {
            self.is_eof = true;
        } else {
            self.buffer.fill(bytes_read);
        }
        Ok(bytes_read)
    }

    fn is_eof(&self) -> bool {
        self.is_eof
    }

    /// Apply a nom parser to the input. Returns Err if the nom parser fails.
    fn run_parser<U>(
        &mut self,
        parser: impl Fn(&[u8]) -> IResult<&[u8], U>,
        detailed_errors: bool,
    ) -> Result<U, StreamParserError> {
        use nom::Context::Code;
        loop {
            let res = match parser(self.buffer.data()) {
                Ok((i, o)) => {
                    let length = self.buffer.data().offset(i);
                    Some((length, o))
                }
                Err(nom::Err::Incomplete(_)) => {
                    // get more data
                    None
                }
                Err(nom::Err::Error(Code(i, e))) | Err(nom::Err::Failure(Code(i, e))) => {
                    let e_slice = if detailed_errors {
                        Some(i[..cmp::min(i.len(), MAX_CONTEXT_BYTES)].to_owned())
                    } else {
                        None
                    };
                    return Err(StreamParserError::StreamParser(e_slice, e));
                }
            };

            match res {
                Some((length, o)) => {
                    self.buffer.consume(length);
                    return Ok(o);
                }
                None => {
                    //refill buffer
                    if self.fill_buffer()? == 0 {
                        return Err(StreamParserError::EOF);
                    }
                }
            }
        }
    }

    /// Try to apply a nom parser, returns None if the parser fails.
    fn try_run_parser<U>(
        &mut self,
        parser: impl Fn(&[u8]) -> IResult<&[u8], U>,
        fail_on_eof: bool,
    ) -> Result<Option<U>, GbParserError> {
        match self.run_parser(parser, false) {
            Ok(o) => Ok(Some(o)),
            Err(StreamParserError::EOF) => {
                if fail_on_eof {
                    Err(StreamParserError::EOF.into())
                } else {
                    Ok(None)
                }
            }
            Err(StreamParserError::StreamParser(_, _)) => Ok(None),
            Err(StreamParserError::Io(e)) => Err(StreamParserError::Io(e).into()),
        }
    }

    /// Apply a nom parser until it fails
    fn run_parser_many0<U>(
        &mut self,
        parser: impl Fn(&[u8]) -> IResult<&[u8], U>,
    ) -> IoResult<Vec<U>> {
        let mut res = Vec::new();
        loop {
            match self.run_parser(&parser, false) {
                Ok(o) => {
                    res.push(o);
                }
                Err(StreamParserError::Io(e)) => {
                    return Err(e);
                }
                _ => {
                    break;
                }
            }
        }
        Ok(res)
    }

    /// Parses the raw sequence data, ignoring whitespace and line numbers
    fn parse_seq_data(&mut self, len: Option<usize>) -> Result<Vec<u8>, GbParserError> {
        let mut s = if let Some(len) = len {
            Vec::with_capacity(cmp::min(len, REASONABLE_SEQ_LEN))
        } else {
            Vec::new()
        };
        loop {
            let mut bytes_read = 0;
            let mut end_of_sequence = false;
            for &b in self.buffer.data() {
                match b {
                    b if b.is_alpha() => {
                        s.push(b);
                    }
                    b'/' => {
                        end_of_sequence = true;
                        break;
                    }
                    b if b.is_dec_digit() => {}
                    b' ' | b'\r' | b'\n' => {}
                    x => {
                        return Err(GbParserError::SyntaxError(format!(
                            "Unexpected char '{}' ({}) in sequence",
                            String::from_utf8_lossy(&[x]), // Only display printable chars
                            x
                        )));
                    }
                }
                bytes_read += 1;
            }
            self.buffer.consume(bytes_read);
            if end_of_sequence {
                // Check that we got everything, if possible
                if let Some(len) = len {
                    if len != s.len() {
                        return Err(GbParserError::SyntaxError(format!(
                            "Got {} bytes of sequence, LOCUS promised {}",
                            s.len(),
                            len
                        )));
                    }
                }
                break;
            }
            if self.fill_buffer()? == 0 {
                if len == Some(s.len()) {
                    warn!("Unexpected EOF while parsing sequence data. Length is correct, continuing.");
                    break;
                } else {
                    // We don't know the length, so we can't know if we have everything
                    return Err(GbParserError::SyntaxError(format!("Unexpected EOF!")));
                }
            }
        }
        Ok(s)
    }

    pub fn read_one_record(&mut self) -> Result<Option<Seq>, GbParserError> {
        // skip preamble such as the header of Genbank .SEQ files
        self.try_run_parser(&skip_preamble, false)?;
        let locus = match self.run_parser(&locus, true) {
            Ok(locus) => locus,
            Err(StreamParserError::EOF) => {
                return Ok(None);
            }
            Err(e) => {
                return Err(e.into());
            }
        };
        let seq = Seq {
            name: locus.name,
            topology: locus.topology,
            date: locus.date,
            len: locus.len,
            molecule_type: locus.molecule_type,
            division: locus.division,
            ..Seq::empty()
        };
        let fields = self.run_parser_many0(&any_field)?;
        let mut seq = fill_seq_fields(seq, fields).map_err(GbParserError::SyntaxError)?; //TODO: Proper error handling
        if self.try_run_parser(&features_header, true)?.is_some() {
            seq.features = self.run_parser_many0(&feature)?;
        }
        self.try_run_parser(&base_count, true)?;
        seq.contig = self.try_run_parser(&contig_text, true)?;
        if self.try_run_parser(&origin_tag, true)?.is_some() {
            seq.seq = self.parse_seq_data(seq.len)?;
        }

        // To be permissive, if we made it this far and it's EOF we'll let the
        // '//' slip
        if self.buffer.empty() && self.is_eof() {
            return Ok(Some(seq));
        }

        self.run_parser(&double_slash, true)?;
        self.run_parser_many0(&line_ending_type_hack)?;
        Ok(Some(seq))
    }
}
