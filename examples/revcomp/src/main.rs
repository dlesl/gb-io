extern crate env_logger;
extern crate log;

#[macro_use]
extern crate clap;
extern crate gb_io;

use std::fs::File;
use std::io;

use gb_io::reader::SeqReader;

fn main() {
    let env = env_logger::Env::default().filter_or(env_logger::DEFAULT_FILTER_ENV, "warn");
    env_logger::Builder::from_env(env).init();
    let args = clap_app!(revcomp =>
        (about: "Reverse complement a genbank file")
        (version: env!("CARGO_PKG_VERSION"))
        (@arg FILE: +required "Genbank file to use as input")
    ).get_matches();
    let file = value_t_or_exit!(args.value_of("FILE"), String);
    let stdout = io::stdout();
    for r in SeqReader::new(File::open(&file).expect("Error reading input file")) {
        r.expect("Error parsing input file")
            .revcomp()
            .write(stdout.lock())
            .expect("Error writing output");
    }
}
