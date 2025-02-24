use std::io;

#[derive(Debug, Error)]
pub enum GbParserError {
    #[error("Syntax error: {0}")]
    SyntaxError(String),
    #[error("{0}")]
    Io(#[source] #[from] io::Error),
}
