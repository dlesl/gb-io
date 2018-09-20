#[macro_use]
extern crate bencher;
extern crate gb_io;


use bencher::Bencher;
use gb_io::reader::{parse_slice, SeqReader};


const ECOLI: &[u8] = include_bytes!("../tests/mg1655.gb");

fn bench_ecoli_slice(b: &mut Bencher) {
    b.iter(|| parse_slice(ECOLI).unwrap());
}

fn bench_ecoli_streaming(b: &mut Bencher) {
    b.iter(|| SeqReader::new(ECOLI).next().unwrap().unwrap());
}

benchmark_group!(benches, bench_ecoli_slice, bench_ecoli_streaming);
benchmark_main!(benches);
