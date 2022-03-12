# gb-io

[![crates.io link](https://img.shields.io/crates/v/gb-io.svg)](https://crates.io/crates/gb-io)
[![docs.rs link](https://docs.rs/gb-io/badge.svg)](https://docs.rs/gb-io/)

This is a library for parsing and working with Genbank (.gb) files written in
Rust. It supports reading, writing, and extracting parts of a sequence while
retaining feature annotations.

It should be able to handle most files out there. Feedback, improvements, and
details of any .gb files it chokes on are welcome!

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

### Python bindings
[Martin Larralde](https://github.com/althonos) has written [Python
bindings](https://pypi.org/project/gb-io/) for `gb-io`'s parser.
