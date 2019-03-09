#[macro_use]
extern crate bencher;
extern crate gb_io;

use bencher::Bencher;
use gb_io::reader::{parse_slice, SeqReader};

const ECOLI: &[u8] = include_bytes!("../tests/mg1655.gb");

fn ecoli_slice(b: &mut Bencher) {
    b.iter(|| parse_slice(ECOLI).unwrap());
}

fn ecoli_streaming(b: &mut Bencher) {
    b.iter(|| SeqReader::new(ECOLI).next().unwrap().unwrap());
}

fn ecoli_revcomp(b: &mut Bencher) {
    let ec = SeqReader::new(ECOLI).next().unwrap().unwrap();
    b.iter(|| ec.revcomp());
}

fn ecoli_set_origin(b: &mut Bencher) {
    let ec = SeqReader::new(ECOLI).next().unwrap().unwrap();
    b.iter(|| ec.set_origin(1_000_000));
}

benchmark_group!(benches, ecoli_slice, ecoli_streaming, ecoli_revcomp, ecoli_set_origin);
benchmark_main!(benches);
