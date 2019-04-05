#[macro_use]
extern crate gb_io;
extern crate bio;

use std::env::args;
use std::fs::File;
use std::str;
use std::io::{BufWriter, Write, stdout};

use gb_io::reader::SeqReader;

fn main() {
    // output in Genbank format?
    let gb_out = args().any(|arg| arg == "--gb");
    let stdout = stdout();
    let mut out = BufWriter::new(stdout.lock());
    for r in SeqReader::new(File::open("tests/mg1655.gb").unwrap()) {
        let r = r.unwrap();
        let record_name = r.name.clone().unwrap();
        for f in &r.features {
            if f.kind == feature_kind!("CDS") {
                let name = format!(
                    "{}_{}_{}",
                    record_name,
                    f.qualifier_values(qualifier_key!("locus_tag"))
                        .next()
                        .unwrap(),
                    f.qualifier_values(qualifier_key!("gene")).next().unwrap()
                );
                if gb_out {
                    let mut extracted = r.extract_location(&f.location).unwrap();
                    extracted.name = Some(name);
                    extracted.write(&mut out).unwrap();
                } else {
                    // FASTA format
                    writeln!(&mut out, ">{}", name).unwrap();
                    let extracted = r.extract_location_seq(&f.location).unwrap();
                    for line in extracted.chunks(60) {
                        writeln!(&mut out, "{}", str::from_utf8(line).unwrap()).unwrap();
                    }
                }
            }
        }
    }
}
