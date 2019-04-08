# gb-io

[![docs.rs link](https://docs.rs/gb-io/badge.svg)](https://docs.rs/gb-io/)

This is a library for working with Genbank (.gb) files written in
Rust. It supports reading, writing, and extracting parts of a sequence
while retaining feature annotations.

The main aim of this project was to learn Rust, and it is not yet
complete, however in its current state it should be able to handle
most files. Feedback, improvements, and details of any .gb files it
chokes on are welcome!

### Example
Reverse complement a sequence, retaining feature annotations.

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
