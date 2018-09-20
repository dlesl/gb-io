use std::io;

#[derive(Debug, Fail)]
pub enum GbParserError {
    #[fail(display = "Syntax error: {}", _0)]
    SyntaxError(String),
    #[fail(display = "{}", _0)]
    Io(#[cause] io::Error),
}

impl From<io::Error> for GbParserError {
    fn from(e: io::Error) -> GbParserError {
        GbParserError::Io(e)
    }
}
