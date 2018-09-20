extern crate flate2;
extern crate gb_io;
extern crate reqwest;

use std::io;
use std::str;

use flate2::read::GzDecoder;
use gb_io::reader::SeqReader;

fn main() {
    let files = (1..10).map(|i| format!("https://ftp.ncbi.nih.gov/genbank/gbmam{}.seq.gz", i));
    for url in files {
        println!("Trying: {}", url);
        let resp = reqwest::get(&url).unwrap();
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
        }
        println!(
            "Found {} records, containing {} bp of which {} were read",
            n_records, total_header_bp, total_seq_bp
        );
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
