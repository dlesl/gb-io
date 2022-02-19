// mapping based on https://github.com/rust-bio/rust-bio/blob/b6cb8699fb7f16e741a7840f5bcc2d850938a37a/src/alphabets/dna.rs
fn comp(base: &u8) -> u8 {
    match base {
        // uppercase
        b'A' => b'T',
        b'G' => b'C',
        b'C' => b'G',
        b'T' => b'A',
        b'Y' => b'R',
        b'R' => b'Y',
        b'W' => b'W',
        b'S' => b'S',
        b'K' => b'M',
        b'M' => b'K',
        b'D' => b'H',
        b'V' => b'B',
        b'H' => b'D',
        b'B' => b'V',
        b'N' => b'N',
        // lowercase
        b'a' => b't',
        b'g' => b'c',
        b'c' => b'g',
        b't' => b'a',
        b'y' => b'r',
        b'r' => b'y',
        b'w' => b'w',
        b's' => b's',
        b'k' => b'm',
        b'm' => b'k',
        b'd' => b'h',
        b'v' => b'b',
        b'h' => b'd',
        b'b' => b'v',
        b'n' => b'n',
        _ => *base,
    }
}

pub fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(comp).collect()
}

#[cfg(test)]
mod test {
    use super::*;
    #[test]
    fn test_revcomp() {
        assert_eq!(b"gAtCnN"[..], revcomp(b"NnGaTc"));
    }
}
