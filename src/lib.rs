//! This is a crate for working with annotated biological sequences stored in
//! "Genbank" format. It focuses on parsing and writing Genbank files and
//! provides some methods for extracting regions from a sequence, preserving
//! annotations and respecting circular molecules.

// #![warn(missing_docs)]
#![cfg_attr(feature = "cargo-clippy", allow(clippy::len_without_is_empty))]
#![cfg_attr(feature = "cargo-clippy", allow(clippy::useless_format))]
#![cfg_attr(feature = "cargo-clippy", allow(clippy::unreadable_literal))]
#![cfg_attr(feature = "cargo-clippy", allow(clippy::upper_case_acronyms))]
#![recursion_limit = "128"]
extern crate circular;
#[macro_use]
extern crate err_derive;
extern crate itertools;
#[macro_use]
extern crate log;
#[macro_use]
extern crate nom;
#[cfg(feature = "serde")]
#[macro_use]
extern crate serde;
extern crate string_cache;

mod errors;

include!(concat!(env!("OUT_DIR"), "/atoms.rs")); // for QualifierKey, FeatureKind

pub mod seq;
pub mod reader;
pub mod writer;
mod dna;

#[cfg(test)]
pub mod tests {
    use crate::errors::GbParserError;
    use crate::reader::*;
    use crate::seq::*;
    use std::fs::File;

    extern crate glob;
    use crate::tests::glob::glob;

    extern crate env_logger;
    pub fn init() {
        let _ = env_logger::Builder::from_default_env().try_init();
    }


    #[test]
    fn streaming() {
        init();
        let f = File::open("tests/mg1655.gb").unwrap();
        let seq = SeqReader::new(f)
            .collect::<Result<Vec<_>, GbParserError>>()
            .unwrap();
        assert_eq!(seq.len(), 1);
        println!("Got {} nts", seq[0].seq.len());
    }
    #[test]
    fn test_seq() {
        init();
        let f = File::open("tests/biopython_tests/gbvrl1_start.seq").unwrap();
        let seqs = SeqReader::new(f);
        for s in seqs {
            let s = s.unwrap();
            println!("{:?}", s.name);
        }
    }
    #[test]
    fn ecoli_read_write() {
        init();
        let ecoli = include_bytes!("../tests/mg1655.gb");
        let r = SeqReader::new(&ecoli[..]).next().unwrap().unwrap();
        let mut out = Vec::new();
        r.write(&mut out).unwrap();
        out.push(b'\n');
        // use std::io::Write;
        // File::create("/tmp/1").unwrap().write_all(&out[..]).unwrap();
        assert_eq!(&ecoli[..], &out[..]);
        assert_eq!(r.len, Some(r.seq.len()));
        assert_eq!(r.features.len(), 9412);
        assert_eq!(r.topology, Topology::Circular);
    }

    #[test]
    fn parse_circular() {
        init();
        let circ = parse_slice(include_bytes!("../tests/circ.gb")).unwrap();
        assert_eq!(circ[0].topology, Topology::Circular);
    }

    #[test]
    fn circular_set_origin() {
        init();
        let circ = parse_slice(include_bytes!("../tests/circ.gb")).unwrap().pop().unwrap();
        for i in 1..10 {
            let rotated = circ.set_origin(i);
            println!("lacZ: {:?}", rotated.features[2].location);
            let rotated_back = rotated.set_origin(circ.len() - i);
            // assert_eq!(rotated_back.seq, circ.seq);
            // assert_eq!(rotated_back.features.len(), circ.features.len());
            // for f in 0..rotated_back.features.len() {
            //     println!("{}", f);
            //     assert_eq!(rotated_back.features[f], circ.features[f]);
            // }
            assert_eq!(rotated_back, circ);
        }
    }
    #[test]

    fn test_multiple_records() {
        init();
        let ls_orchid =
            parse_slice(include_bytes!("../tests/biopython_tests/ls_orchid.gb")).unwrap();
        assert_eq!(ls_orchid.len(), 94);
    }

    #[test]
    fn biopython_tests() {
        init();
        for f in glob("tests/biopython_tests/*.gb").unwrap() {
            let f = f.unwrap();
            println!("Testing: {:?}", f);
            let records = SeqReader::new(::std::fs::File::open(f).unwrap());
            for r in records {
                let r = r.unwrap();
                println!("{:?}", r.name);
                if let Some(len) = r.len {
                    if r.contig.is_none() || !r.seq.is_empty() {
                        assert_eq!(len, r.seq.len());
                    }
                }
            }
        }
    }
}
