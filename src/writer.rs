use itertools::Itertools;
use crate::seq::{Date, QualifierKey, Seq};
use std::convert::AsRef;
use std::io::{self, Write};

// Following constants taken from Biopython's InsdcIO.py
const MAX_WIDTH: usize = 79;
const QUALIFIER_INDENT: &str = "                     ";
const FIELD_INDENT: &str = "            ";
// Commented out some of these, need to check if they're still in
// use
const FTQUAL_NO_QUOTE: &[QualifierKey] = &[
    qualifier_key!("anticodon"),
    //qualifier_key!("citation"),
    qualifier_key!("codon_start"),
    //qualifier_key!("compare"),
    //qualifier_key!("direction"),
    qualifier_key!("estimated_length"),
    //qualifier_key!("mod_base"),
    qualifier_key!("number"),
    qualifier_key!("rpt_type"),
    //qualifier_key!("rpt_unit_range"),
    //qualifier_key!("tag_peptide"),
    qualifier_key!("transl_except"),
    qualifier_key!("transl_table"),
];

const POS_QUAL: &[QualifierKey] = &[
    // should be formatted like Locations
    qualifier_key!("transl_except"),
    qualifier_key!("anticodon"),
];



#[derive(Debug)]
pub struct SeqWriter<W: Write> {
    stream: W,
    truncate_locus: bool,
    escape_locus: bool,
}

impl<W: Write> SeqWriter<W> {
    /// Create a new `SeqWriter` to write GenBank files to the given stream.
    pub fn new(stream: W) -> Self {
        Self {
            stream,
            truncate_locus: false,
            escape_locus: true,
        }
    }

    /// Set the behaviour regarding locus name truncation.
    ///
    /// If `true` (the default when creating a `SeqWriter`, for backwards
    /// compatibility), the locus fields (`name` and `molecule_type`) will be
    /// truncated if they are too long so that the LOCUS line is no longer
    /// than 79 characters.
    ///
    /// By setting `truncate` to `false`, the full strings are written to the
    /// LOCUS line, resulting in fields that may not in usual positions if
    /// the name or molecule types are too long. Most parsers (including `gb-io`
    /// and Biopython) should however be able to process these.
    pub fn truncate_locus(&mut self, truncate: bool) -> &mut Self {
        self.truncate_locus = truncate;
        self
    }

    /// Set the behaviour regarding locus name escaping.
    ///
    /// If `true` (the default when creating a `SeqWriter`, for backwards
    /// compatibility), any whitespace in the locus will be escaped with
    /// and underscore (_).
    ///
    /// Set to `false` to write the locus as-is, which could result in an
    /// invalid LOCUS line that some parsers may not be able to process.
    pub fn escape_locus(&mut self, escape: bool) -> &mut Self {
        self.escape_locus = escape;
        self
    }

    /// Generate the locus line for the record.
    ///
    /// Ported from Biopython (InsdcIO.py).
    fn locus_line(&self, record: &Seq) -> String {
        let mut locus = record.name.clone().unwrap_or_else(|| {
            record
                .accession
                .clone()
                .unwrap_or_else(|| "UNTITLED".into())
        });
        let length = format!("{}", record.len());
        if self.truncate_locus && locus.len() + 1 + length.len() > 28 {
            locus = locus[..27 - length.len()].into();
        }

        if self.escape_locus && locus.split_whitespace().count() > 1 {
            locus = locus.split_whitespace().join("_");
        }

        let units = "bp"; // TODO: 'aa' support

        let mol_type = match record.molecule_type {
            Some(ref m) if m.len() > 7 && self.truncate_locus => "",
            Some(ref m) => m,
            None => ""
        };

        if locus.len() + 1 + length.len() >= 28 {
            let length = format!(" {}", length);
            locus.push_str(&length);
        } else {
            let length = format!("{:>28}", length);
            let rest = &length[locus.len()..];
            locus.push_str(rest);
        }

        format!(
            "LOCUS       {} {}    {:<7} {:<8} {} {}\n",
            locus,
            units,
            mol_type,
            record.topology,
            record.division,
            record
                .date
                .as_ref()
                .unwrap_or(&Date::from_ymd(1970, 1, 1).unwrap())
        )
    }

    /// Write the sequence to the stream.
    pub fn write(&mut self, record: &Seq) -> io::Result<()> {
        // LOCUS

        let locus_line = self.locus_line(record);
        write!(&mut self.stream, "{}", locus_line)?;

        // Fields

        write_field_maybe(&mut self.stream, &record.definition, "DEFINITION")?;
        write_field_maybe(&mut self.stream, &record.accession, "ACCESSION")?;
        write_field_maybe(&mut self.stream, &record.version, "VERSION")?;
        write_field_maybe(&mut self.stream, &record.dblink, "DBLINK")?;
        write_field_maybe(&mut self.stream, &record.keywords, "KEYWORDS")?;
        if let Some(ref source) = record.source {
            write_field(&mut self.stream, &source.source, "SOURCE")?;
            write_field_maybe(&mut self.stream, &source.organism, "  ORGANISM")?;
        }
        for r in &record.references {
            write_field(&mut self.stream, &r.description, "REFERENCE")?;
            for a in &r.authors {
                write_field(&mut self.stream, a, "  AUTHORS")?;
            }
            write_field_maybe(&mut self.stream, &r.consortium, "  CONSRTM")?;
            write_field(&mut self.stream, &r.title, "  TITLE")?;
            write_field_maybe(&mut self.stream, &r.journal, "  JOURNAL")?;
            write_field_maybe(&mut self.stream, &r.pubmed, "   PUBMED")?; // TODO: this should be a subfield of Journal
            write_field_maybe(&mut self.stream, &r.remark, "  REMARK")?;
        }
        for comment in &record.comments {
            write_field(&mut self.stream, comment, "COMMENT")?;
        }

        // Features

        if !record.features.is_empty() {
            self.stream.write_all(b"FEATURES             Location/Qualifiers\n")?;
            for f in &record.features {
                let first_indent = format!("     {:<15} ", f.kind);
                let location = f.location.to_gb_format();
                wrap_location(
                    &mut self.stream,
                    &location,
                    MAX_WIDTH,
                    first_indent.as_str(),
                    QUALIFIER_INDENT,
                )?;
                for &(ref key, ref val) in &f.qualifiers {
                    match *val {
                        None => writeln!(&mut self.stream, "{}/{}", QUALIFIER_INDENT, key)?,
                        Some(ref val) => {
                            let quote = !FTQUAL_NO_QUOTE.iter().any(|x| x == key);
                            let first_indent = format!("{}/{}=", QUALIFIER_INDENT, key);
                            if POS_QUAL.iter().any(|x| x == key) {
                                wrap_location(
                                    &mut self.stream,
                                    val,
                                    MAX_WIDTH,
                                    first_indent.as_str(),
                                    QUALIFIER_INDENT,
                                )?;
                            } else {
                                wrap_text(
                                    &mut self.stream,
                                    val,
                                    MAX_WIDTH,
                                    first_indent.as_str(),
                                    QUALIFIER_INDENT,
                                    quote,
                                )?;
                            }
                        }
                    }
                }
            }
        }

        // CONTIG, maybe

        if let Some(ref contig) = record.contig {
            wrap_location(
                &mut self.stream,
                &contig.to_gb_format(),
                MAX_WIDTH,
                "CONTIG      ",
                FIELD_INDENT,
            )?;
        }

        // ORIGIN

        if !record.seq.is_empty() {
            let mut line = Vec::with_capacity(79);
            write!(&mut line, "ORIGIN      ")?;
            for (i, &b) in record.seq.iter().enumerate() {
                if i % 60 == 0 {
                    line.push(b'\n');
                    self.stream.write_all(&line)?;
                    line.clear();
                    write!(&mut line, "{:>9}", i + 1)?;
                }
                if i % 10 == 0 {
                    line.push(b' ');
                }
                line.push(b);
            }
            line.push(b'\n');
            self.stream.write_all(&line)?;
        }
        writeln!(&mut self.stream, "//")?;
        Ok(())
    }
}

pub fn write<T: Write>(file: T, record: &Seq) -> io::Result<()> {
    SeqWriter::new(file).write(record)
}

fn write_field<T>(mut file: T, field: &str, keyword: &str) -> io::Result<()>
where
    T: Write,
{
    let keyword = format!("{:<12}", keyword);
    wrap_text(
        &mut file,
        field,
        MAX_WIDTH,
        keyword.as_str(),
        FIELD_INDENT,
        false,
    )
}

fn write_field_maybe<T, U>(mut file: T, field: &Option<U>, keyword: &str) -> io::Result<()>
where
    T: Write,
    U: AsRef<str>,
{
    if let Some(ref field) = *field {
        write_field(&mut file, field.as_ref(), keyword)?;
    };
    Ok(())
}

/// Wrap a genbank location, splitting on commas if possible
fn wrap_location<T: Write>(
    mut file: T,
    mut text: &str,
    max_width: usize,
    first_indent: &str,
    subsequent_ident: &str,
) -> io::Result<()> {
    let mut indent = first_indent;
    while indent.len() + text.len() > max_width {
        let split_at = if let Some(location) = text[..max_width - indent.len()].rfind(',') {
            location + 1
        } else {
            warn!("Couldn't split location appropriately!");
            max_width - indent.len()
        };
        writeln!(file, "{}{}", indent, &text[..split_at])?;
        text = &text[split_at..];
        indent = subsequent_ident;
    }
    if !text.is_empty() {
        writeln!(file, "{}{}", indent, text)?;
    }
    Ok(())
}

/// Wrap a line of text in a Genbank file. The line is appended to `line`, which
/// may already contain up to `max_len` - 1 bytes. Wrapping occurs on a space if
/// possible. Outputs valid UTF-8, however lines are wrapped to their length in
/// bytes, not characters, and UTF-8 characters that consist of multiple `char`s
/// are likely to get mangled. Returns the rest of `input`. If `quote` is set,
/// `"` will be escaped with `""`
fn wrap_get_line<'a>(line: &mut String, input: &'a str, max_len: usize, quote: bool) -> &'a str {
    assert!(line.len() < max_len);
    let mut consumed = 0;
    let mut last_space_in = None;
    let mut last_space_out = 0; // need to keep track of this separately because
                                // of escape sequences
    let mut i = input.char_indices();
    while line.len() < max_len {
        if let Some((idx, ch)) = i.next() {
            match ch {
                ' ' => {
                    last_space_in = Some(idx);
                    last_space_out = line.len();
                    line.push(' ');
                }
                '\n' => {
                    return &input[idx + 1..];
                }
                '"' if quote => {
                    // if we're one character away from the end of the line,
                    // wrap now to avoid splitting the escape sequence
                    if line.len() >= max_len - 1 {
                        break;
                    } else {
                        line.push_str("\"\"");
                    }
                }
                _ => {
                    line.push(ch);
                }
            }
            consumed = idx + ch.len_utf8();
        } else {
            // end of input
            return &input[consumed..];
        }
    }
    // So it's time to split the line. First check if it's really necessary,
    // maybe the next character is a newline or space or the end
    if let Some((idx_next, next)) = i.next() {
        if next == ' ' || next == '\n' {
            return &input[idx_next + 1..];
        }
    } else {
        // it's the end, no need to wrap
        assert!(consumed == input.len());
        return &input[consumed..];
    }
    // try to wrap at last space
    if let Some(last_space_in) = last_space_in {
        line.truncate(last_space_out);
        &input[last_space_in + 1..]
    } else {
        // there's no way to split the line nicely
        &input[consumed..]
    }
}

fn wrap_text<T: Write>(
    mut file: T,
    mut text: &str,
    max_width: usize,
    first_indent: &str,
    subsequent_indent: &str,
    quote: bool,
) -> io::Result<()> {
    let mut line = String::with_capacity(max_width);
    let mut indent = first_indent;
    if quote {
        line.push('"');
    }
    text = wrap_get_line(&mut line, text, max_width - indent.len(), quote);
    write!(file, "{}{}", indent, line)?;
    while !text.is_empty() {
        indent = subsequent_indent;
        line.clear();
        text = wrap_get_line(&mut line, text, max_width - subsequent_indent.len(), quote);
        write!(file, "\n{}{}", subsequent_indent, line)?;
    }
    if quote {
        if line.len() + indent.len() >= max_width {
            writeln!(file, "\n{}\"", subsequent_indent)?;
        } else {
            writeln!(file, "\"")?;
        }
    } else {
        writeln!(file)?;
    }
    Ok(())
}

#[cfg(test)]
pub mod tests {

    use super::*;
    use std::io::BufRead;

    #[test]
    fn truncate_locus() {
        let mut seq = Seq::empty();
        seq.name = Some(String::from("1122217.SAMN02441331.KB899611_14"));

        let mut out = Vec::new();
        SeqWriter::new(&mut out).truncate_locus(true).write(&seq).unwrap();

        let mut line = String::new();
        std::io::Cursor::new(&out).read_line(&mut line).unwrap();
        assert_eq!(
            line,
            "LOCUS       1122217.SAMN02441331.KB899 0 bp            linear UNK 01-JAN-1970\n"
        );
    }

    #[test]
    fn no_truncate_locus() {
        let mut seq = Seq::empty();
        seq.name = Some(String::from("1122217.SAMN02441331.KB899611_14"));

        let mut out = Vec::new();
        SeqWriter::new(&mut out).truncate_locus(false).write(&seq).unwrap();

        let mut line = String::new();
        std::io::Cursor::new(&out).read_line(&mut line).unwrap();
        assert_eq!(
            line,
            "LOCUS       1122217.SAMN02441331.KB899611_14 0 bp            linear UNK 01-JAN-1970\n"
        );
    }
}
