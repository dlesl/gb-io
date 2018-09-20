# gb-io
This is a library for working with Genbank (.gb) files written in
Rust. It aims to support reading, writing, and extracting parts of a
sequence while retaining feature annotations.

The main aim of this project was to learn Rust, and the code is far
from complete (or nice), however it should work for the most
part. Feedback, improvements, and details of any .gb files it chokes
on are welcome!

### Example
```rust
extern crate gb_io;

use std::fs::File;
use std::io;

use gb_io::reader::SeqReader;

fn main() {
    let file = File::open("mg1655.gb").unwrap();
    let stdout = io::stdout();
    for seq in SeqReader::new(file) {
        let seq = seq.unwrap();
        let rc = seq.revcomp();
        rc.write(stdout.lock()).unwrap();
    }
}
```
