extern crate flate2;
extern crate gb_io;
extern crate reqwest;
extern crate itertools;

use std::borrow::Cow;
use std::collections::HashMap;
use std::{io, thread};
use std::str;
use std::sync::{Arc, Mutex};
use flate2::read::GzDecoder;
use itertools::Itertools;
use gb_io::reader::SeqReader;

fn main() {
    let mut files = Vec::new();
    for i in 1..10 {
        files.push(format!("https://ftp.ncbi.nlm.nih.gov/genbank/gbbct{}.seq.gz", i));
        files.push(format!("https://ftp.ncbi.nlm.nih.gov/genbank/gbmam{}.seq.gz", i));
    }
    let stats = Arc::new(Mutex::new(Stats::default()));
    let mut threads = Vec::new();
    for url in files {
        let stats = stats.clone();
        threads.push(thread::spawn(move || {
            println!("Trying: {}", url);
            let resp = reqwest::blocking::get(&url).unwrap();
            assert!(resp.status().is_success());
            let data = GzDecoder::new(resp);
            let mut total_seq_bp = 0;
            let mut total_header_bp = 0;
            let mut n_records = 0;
            for s in SeqReader::new(Latin1ToUtf8(data)) {
                let s = s.unwrap();
                let header_len = s.len.unwrap_or(0);
                let seq_len = s.seq.len();
                total_header_bp += header_len;
                total_seq_bp += seq_len;
                n_records += 1;
                let mut stats_mut = stats.lock().unwrap();
                for f in s.features {
                    *stats_mut.feature_kind_frequencies.entry(f.kind).or_insert(0) += 1;
                    for (q, _) in f.qualifiers {
                        *stats_mut.qualifier_key_frequencies.entry(q).or_insert(0) += 1;
                    }
                }
            }
            println!(
                "{}: Found {} records, containing {} bp of which {} were read",
                url, n_records, total_header_bp, total_seq_bp
            );
            assert_eq!(total_header_bp, total_seq_bp);
            let stats_mut = stats.lock().unwrap();
            println!("Top feature kinds so far:");
            print_percentile(&stats_mut.feature_kind_frequencies, ' ', 0.99);
            println!("Top qualifier kinds so far:");
            print_percentile(&stats_mut.qualifier_key_frequencies, '=', 0.99);
        }));
    }
    for t in threads {
        t.join().unwrap();
    }
}

#[derive(Default)]
struct Stats {
    qualifier_key_frequencies: HashMap<Cow<'static, str>, usize>,
    feature_kind_frequencies: HashMap<Cow<'static, str>, usize>
}

// used to generate the match statements in feature_table.rs
fn print_percentile(m: &HashMap<Cow<str>, usize>, next_char: char, percentile: f64) {
    let total_count: usize = m.iter().map(|(_, v)| v).sum();
    let count_to_take = (total_count as f64 * percentile) as usize;
    let mut taken = 0;
    for (k, v) in m.iter().sorted_by_key(|(_, v)| *v).rev() {
        println!(
            "[{}, b'{}', ..] => Ok((&input[{}..], Cow::from(\"{}\"))),",
            k.chars().map(|c| format!("b'{}'", c)).join(", "),
            next_char,
            k.chars().count(),
            k
        );
        taken += v;
        if taken > count_to_take {
            return;
        }
    }
}

// Some of the records in GenBank contain invalid latin1-encoded characters, for
// example '®' after the name of the software used to generate the sequence.
// Since we're strict about only allowing UTF-8, we need to convert these
// symbols before parsing them. We also log the problem, which should be
// reported to NCBI, see https://github.com/biopython/biopython/issues/1361

struct Latin1ToUtf8<T: io::Read>(T);

impl<T: io::Read> io::Read for Latin1ToUtf8<T> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        // This is not very efficient, but the result of the decoding could be
        // up to 2x the length of the original data (if every byte is invalid).

        // Since these characters are very rare, we first read into the provided
        // buffer. If it contains an invalid character, we allocate a new buffer
        // to do the expansion

        let max_len = buf.len() / 2;
        let res = self.0.read(&mut buf[..max_len])?;

        match str::from_utf8(&buf[..res]) {
            Ok(_) => Ok(res),
            Err(_) => {
                println!("*************** ENCODING ERROR! ******************");
                let new_buf: String = buf[..res].iter().cloned().map(char::from).collect();
                let bytes = new_buf.as_bytes();
                buf[..bytes.len()].copy_from_slice(&bytes);
                Ok(bytes.len())
            }
        }
    }
}

#[test]
fn latin1_to_utf8() {
    use std::io::Read;
    let ascii = b"hello";
    let mut decoded = String::new();
    Latin1ToUtf8(&ascii[..])
        .read_to_string(&mut decoded)
        .unwrap();
    assert_eq!(&decoded, "hello");
    let latin1 = b"tjabba tjena hall\xE5";
    decoded = String::new();
    Latin1ToUtf8(&latin1[..])
        .read_to_string(&mut decoded)
        .unwrap();
    assert_eq!(&decoded, "tjabba tjena hallå");
}
