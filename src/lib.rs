//! This is a crate for working with annotated biological sequences stored in
//! "Genbank" format. It focuses on parsing and writing Genbank files and
//! provides some methods for extracting regions from a sequence, preserving
//! annotations and respecting circular molecules.

// #![warn(missing_docs)]
#![cfg_attr(feature = "cargo-clippy", feature(tool_lints))]
#![cfg_attr(feature = "cargo-clippy", allow(clippy::len_without_is_empty))]
#![cfg_attr(feature = "cargo-clippy", allow(clippy::useless_format))]
#![cfg_attr(feature = "cargo-clippy", allow(clippy::unreadable_literal))]
#![recursion_limit = "128"]
extern crate circular;
#[macro_use]
extern crate failure;
extern crate itertools;
#[macro_use]
extern crate log;
#[macro_use]
extern crate nom;
#[cfg(feature = "serde")]
#[macro_use]
extern crate serde;
extern crate string_cache;
extern crate bio;

mod errors;

include!(concat!(env!("OUT_DIR"), "/atoms.rs")); // for QualifierKey, FeatureKind

pub mod seq;
pub mod reader;
mod writer;

#[cfg(test)]
mod tests {
    use errors::GbParserError;
    use reader::*;
    use seq::*;
    use std::fs::File;

    extern crate glob;
    use tests::glob::glob;


    #[test]
    fn streaming() {
        let f = File::open("tests/mg1655.gb").unwrap();
        let seq = SeqReader::new(f)
            .collect::<Result<Vec<_>, GbParserError>>()
            .unwrap();
        assert_eq!(seq.len(), 1);
        println!("Got {} nts", seq[0].seq.len());
    }
    #[test]
    fn test_seq() {
        let f = File::open("tests/biopython_tests/gbvrl1_start.seq").unwrap();
        let seqs = SeqReader::new(f);
        for s in seqs {
            let s = s.unwrap();
            println!("{:?}", s.name);
        }
    }
    #[test]
    fn ecoli_read_write() {
        let ecoli = include_bytes!("../tests/mg1655.gb");
        let r = SeqReader::new(&ecoli[..]).next().unwrap().unwrap();
        let mut out = Vec::new();
        r.write(&mut out).unwrap();
        out.push(b'\n');
        // File::create("/tmp/1").unwrap().write_all(&out[..]).unwrap();
        assert_eq!(&ecoli[..], &out[..]);
        assert_eq!(r.len, Some(r.seq.len()));
        assert_eq!(r.features.len(), 9412);
        assert_eq!(r.topology, Topology::Circular);
    }

    #[test]
    fn parse_circular() {
        let circ = parse_slice(include_bytes!("../tests/circ.gb")).unwrap();
        assert_eq!(circ[0].topology, Topology::Circular);
    }
    #[test]

    fn test_multiple_records() {
        let ls_orchid =
            parse_slice(include_bytes!("../tests/biopython_tests/ls_orchid.gb")).unwrap();
        assert_eq!(ls_orchid.len(), 94);
    }

    #[test]
    fn biopython_tests() {
        for f in glob("tests/biopython_tests/*.gb").unwrap() {
            // TODO: Add the other extensions
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
