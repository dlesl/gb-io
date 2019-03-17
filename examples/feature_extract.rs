#[macro_use]
extern crate gb_io;
extern crate bio;

use std::fs::File;
use std::str;

use gb_io::reader::SeqReader;

fn main() {
    for r in SeqReader::new(File::open("tests/mg1655.gb").unwrap()) {
        let r = r.unwrap();
        let record_name = r.name.clone().unwrap();
        for f in &r.features {
            if f.kind == feature_kind!("CDS") {
                println!(
                    ">{}:{}:{}",
                    record_name,
                    f.qualifier_values(qualifier_key!("locus_tag")).next().unwrap(),
                    f.qualifier_values(qualifier_key!("gene")).next().unwrap()
                );
                let extracted = r.extract_location_seq(&f.location).unwrap();
                for line in extracted.chunks(60) {
                    println!("{}", str::from_utf8(line).unwrap());
                }
            }
        }
    }
}
