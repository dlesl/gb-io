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

// Copied from nom and modified
macro_rules! map_res_custom_error (
  // Internal parser, do not use directly
  (__impl $i:expr, $e:expr, $submac:ident!( $($args:tt)* ), $submac2:ident!( $($args2:tt)* )) => (
    {
      use ::std::result::Result::*;
      use ::nom::Err;

      let i_ = $i.clone();
      ($submac!(i_, $($args)*)).and_then(|(i,o)| {
        match $submac2!(o, $($args2)*) {
          Ok(output) => Ok((i, output)),
          Err(_) => {
            let e = ::nom::ErrorKind::Custom(u32::from($e));
            Err(Err::Error(error_position!($i, e)))
          },
        }
      })
    }
  );
  ($i:expr, $e:expr, $submac:ident!( $($args:tt)* ), $g:expr) => (
    map_res_custom_error!(__impl $i, $e, $submac!($($args)*), call!($g))
  );
  ($i:expr, $e:expr, $submac:ident!( $($args:tt)* ), $submac2:ident!( $($args2:tt)* )) => (
    map_res_custom_error!(__impl $i, $e, $submac!($($args)*), $submac2!($($args2)*))
  );
  ($i:expr, $e:expr, $f:expr, $g:expr) => (
    map_res_custom_error!(__impl $i, $e, call!($f), call!($g))
  );
  ($i:expr, $e:expr, $f:expr, $submac:ident!( $($args:tt)* )) => (
    map_res_custom_error!(__impl $i, call!($f), $submac!($($args)*))
  )
);

macro_rules! to_str (
    ($i:expr, $submac:ident!( $($args:tt)* )) => (
        map_res_custom_error!($i, NomParserError::Utf8Error,
                      $submac!($($args)*),
                               str::from_utf8
        )
    );
    ($i:expr, $f:expr) => (
        to_str!($i, call!($f))
    )
);

macro_rules! to_string (
    ($i:expr, $submac:ident!( $($args:tt)* )) => (
        map_res_custom_error!($i, NomParserError::Utf8Error,
                            $submac!($($args)*),
                            String::from_utf8
        )
    );
    ($i:expr, $f:expr) => (
        to_str!($i, call!($f))
    )
);
