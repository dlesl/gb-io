#[macro_use]
extern crate gb_io;
extern crate bio;

use std::fs::File;
use std::str;

use bio::alphabets::dna::revcomp;
use gb_io::reader::SeqReader;
use gb_io::seq::Position;

fn main() {
    for r in SeqReader::new(File::open("tests/mg1655.gb").unwrap()) {
        let r = r.unwrap();
        let record_name = r.name.clone().unwrap();
        for f in &r.features {
            if f.kind == feature_kind!("CDS") {
                println!(
                    ">{}:{}:{}",
                    record_name,
                    f.get_qualifier_values(&qualifier_key!("locus_tag"))[0],
                    f.get_qualifier_values(&qualifier_key!("gene"))[0]
                );
                let (start, end) = f.pos.find_bounds().unwrap();
                let extracted = match f.pos {
                    Position::Complement(_) => {
                        revcomp(&r.extract_range_seq(start, end + 1)[..]).into()
                    }
                    _ => r.extract_range_seq(start, end + 1),
                };
                for line in extracted.chunks(60) {
                    println!("{}", str::from_utf8(line).unwrap());
                }
            }
        }
    }
}
