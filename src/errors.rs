use std::io;

#[derive(Debug, Error)]
pub enum GbParserError {
    #[error(display = "Syntax error: {}", _0)]
    SyntaxError(String),
    #[error(display = "{}", _0)]
    Io(#[cause] io::Error),
}
