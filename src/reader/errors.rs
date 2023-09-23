/// Error types and some helper macros for using them in nom parsers
use std::fmt;

// At the moment it's very difficult to use custom error types in nom, since
// macros like named_args! don't support it, so we'll stick to u32 until that
// work is completed

// TODO: Messages that explain what feature of the file was unable to be parsed
#[derive(Debug)]
pub enum NomParserError {
    Utf8Error,
    Date,
    Location,
    SequenceLength,
    Unknown,
}

impl From<NomParserError> for u32 {
    fn from(i: NomParserError) -> u32 {
        i as u32
    }
}

impl From<u32> for NomParserError {
    fn from(i: u32) -> NomParserError {
        match i {
            0 => NomParserError::Utf8Error,
            1 => NomParserError::Date,
            2 => NomParserError::Location,
            3 => NomParserError::SequenceLength,
            _ => NomParserError::Unknown,
        }
    }
}

impl fmt::Display for NomParserError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let msg = match *self {
            NomParserError::Utf8Error => "Invalid UTF-8 encountered",
            NomParserError::Date => "Invalid date format",
            NomParserError::Location => "Failed parsing location specifier",
            NomParserError::SequenceLength => {
                "Sequence length didn't match value specified by LOCUS line"
            }
            NomParserError::Unknown => "Unknown parse error", // This should make people happy!
        };
        write!(f, "{}", msg)
    }
}
