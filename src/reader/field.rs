use std::iter::once;
use std::str;

use itertools::intersperse;
use nom::branch::alt;
use nom::bytes::streaming::{is_a, tag};
use nom::character::streaming::{line_ending, not_line_ending, char};
use nom::combinator::{map, map_res, opt, value};
use nom::multi::many0;
use nom::sequence::{delimited, terminated};
use nom::{IResult, Parser};
use nom::multi::fold;

use crate::reader::errors::NomParserError;
use crate::reader::location::location;
use crate::reader::misc;
use crate::seq::{Location, Reference, Seq, Source};

pub fn fields(i: &[u8]) -> IResult<&[u8], Vec<Field>> {
    many0(field).parse(i)
}

/// Matches `indent` spaces
pub fn space_indent<'a>(indent: usize) -> impl Parser<&'a [u8], Output = (), Error = nom::error::Error<&'a [u8]>> {   //-> impl FnMut(&[u8]) -> IResult<&[u8], ()> {
    fold(indent, char(' '), || (), |_, _| () )
}

pub fn field_string<'a>(
    indent: usize,
    name: &'static str,
    keep_ws: bool,
) -> impl Parser<&'a [u8], Output = String, Error = nom::error::Error<&'a [u8]>> {
    map_res(field_bytes(indent, name, keep_ws), String::from_utf8)
}

/// Concatenates lines, optionally interpolating with newlines
fn concat_lines<'a, T: Iterator<Item = &'a [u8]>>(lines: T, keep_linebreaks: bool) -> Vec<u8> {
    if keep_linebreaks {
        intersperse(lines, b"\n").flatten().cloned().collect()
    } else {
        lines.flatten().cloned().collect()
    }
}

// Matches an entry beginning with an (optionally) indented identifier followed
// by one of more indented lines of free text
fn field_bytes(
    indent: usize,
    name: &'static str,
    keep_ws: bool,
) -> impl FnMut(&[u8]) -> IResult<&[u8], Vec<u8>> {
    move |i| {
        let (i, _) = space_indent(indent).parse(i)?;
        let (i, _) = tag(name)(i)?;
        let (i, spaces) = map(is_a(" "), <[_]>::len).parse(i)?;
        let (i, first_line) = terminated(not_line_ending, line_ending).parse(i)?;
        let (i, more_lines) = many0(delimited(
            space_indent(indent + spaces + name.len()),
            not_line_ending,
            line_ending,
        )).parse(i)?;
        Ok((
            i,
            concat_lines(once(first_line).chain(more_lines), keep_ws),
        ))
    }
}

// I think this is deprecated but sometimes it comes after features,
// let's just ignore it.
pub fn base_count(i: &[u8]) -> IResult<&[u8], ()> {
    value((), field_string(0, "BASE COUNT", false)).parse(i)
}

// REFERENCE

fn reference(i: &[u8]) -> IResult<&[u8], Reference> {
    let (i, description) = field_string(0, "REFERENCE", true).parse(i)?;
    let (i, authors) = opt(field_string(2, "AUTHORS", true)).parse(i)?;
    let (i, consortium) = opt(field_string(2, "CONSRTM", true)).parse(i)?;
    let (i, title) = field_string(2, "TITLE", true).parse(i)?;
    let (i, journal) = opt(field_string(2, "JOURNAL", true)).parse(i)?;
    let (i, pubmed) = opt(field_string(3, "PUBMED", false)).parse(i)?;
    let (i, remark) = opt(field_string(2, "REMARK", true)).parse(i)?;
    Ok((
        i,
        Reference {
            description,
            authors,
            consortium,
            title,
            journal,
            pubmed,
            remark,
        },
    ))
}

// These are in CAPS so we can stringify! them and match the tags directly
#[derive(Debug)]
pub enum Field {
    DEFINITION(String),
    ACCESSION(String),
    VERSION(String),
    DBLINK(String),
    KEYWORDS(String),
    SOURCE(Source),
    REFERENCE(Reference),
    COMMENT(String),
    UnrecognisedLine(String),
}

macro_rules! parse_field {
    ($field: ident, $keep_ws:expr) => {
        map(field_string(0, stringify!($field), $keep_ws), Field::$field)
    };
}

fn source(i: &[u8]) -> IResult<&[u8], Field> {
    let (i, source) = field_string(0, "SOURCE", true).parse(i)?;
    let (i, organism) = opt(field_string(2, "ORGANISM", true)).parse(i)?;
    Ok((i, Field::SOURCE(Source { source, organism })))
}

pub fn field(i: &[u8]) -> IResult<&[u8], Field> {
    alt((
        parse_field!(DEFINITION, true),
        parse_field!(ACCESSION, true),
        parse_field!(VERSION, true),
        parse_field!(DBLINK, true),
        parse_field!(KEYWORDS, true),
        source,
        map(reference, Field::REFERENCE),
        parse_field!(COMMENT, true),
        //TODO: unrecognised lines?
        map(misc::ignored_line, |line| {
            Field::UnrecognisedLine(line.into())
        }),
    )).parse(i)
}

/// Create a new Seq with metadata from a Vec<Field>
pub fn fill_seq_fields(mut seq: Seq, fields: Vec<Field>) -> Result<Seq, String> {
    // TODO: Use real errors once we have a way to return them through nom
    // helper function
    fn set_once<T>(field: &str, val: &mut Option<T>, newval: T) -> Result<(), String> {
        match *val {
            Some(_) => Err(format!("Field '{}' occurred twice!", field)),
            None => {
                *val = Some(newval);
                Ok(())
            }
        }
    }

    for item in fields {
        match item {
            Field::DEFINITION(item) => set_once("DEFINITION", &mut seq.definition, item)?,
            Field::ACCESSION(item) => {
                seq.accession = Some(item);
            }
            Field::VERSION(item) => {
                seq.version = Some(item);
            }
            Field::DBLINK(item) => {
                seq.dblink = Some(item);
            }
            Field::KEYWORDS(item) => {
                seq.keywords = Some(item);
            }
            Field::SOURCE(item) => {
                seq.source = Some(item);
            }
            Field::REFERENCE(item) => {
                seq.references.push(item);
            }
            Field::COMMENT(item) => {
                seq.comments.push(item);
            }
            Field::UnrecognisedLine(line) => {
                warn!("Unable to parse: {}", line);
            }
        }
    }
    Ok(seq)
}

pub fn contig_text(i: &[u8]) -> IResult<&[u8], Location> {
    map_res(field_bytes(0, "CONTIG", false), |bytes| {
        match location(&bytes) {
            Ok((_, l)) => Ok(l),
            Err(_) => Err(NomParserError::Location),
        }
    }).parse(i) // TODO: proper error reporting
}

#[cfg(test)]
mod test {
    use super::*;

    fn incomplete_test<T: std::fmt::Debug>(
        input: &[u8],
        remainder: &[u8],
        parser: impl Fn(&[u8]) -> IResult<&[u8], T>,
    ) {
        match parser(input) {
            Ok((i, _)) => {
                assert_eq!(i, remainder);
            }
            x => panic!("{:#?}", x),
        }
        for n in 0..(input.len() - remainder.len()) {
            match parser(&input[..n]) {
                Err(nom::Err::Incomplete(_)) => {}
                x => panic!("{:?} => {:#?}", &input[..n], x),
            }
        }
    }

    #[test]
    fn test_field() {
        let f = br#"COMMENT     REFSEQ INFORMATION: The reference sequence is identical to
            GL000055.2.
            On or before Feb 3, 2014 this sequence version replaced
            gi:224514773, gi:338858164, gi:338858163, gi:224514859.
            Assembly Name: GRCh38.p7 Primary Assembly
            The DNA sequence is composed of genomic sequence, primarily
            finished clones that were sequenced as part of the Human Genome
            Project. PCR products and WGS shotgun sequence have been added
            where necessary to fill gaps or correct errors. All such additions
            are manually curated by GRC staff. For more information see:
            http://genomereference.org.
            This sequence is part of a chromosome or linkage group.  Current
            annotation for this sequence is available on the chromosome or
            linkage group record.
FEATURES             Location/Qualifiers"#;
        incomplete_test(f, b"FEATURES             Location/Qualifiers", field);
    }

    #[test]
    fn chr1_contig() {
        let s = r#"CONTIG      join(gap(10000),NT_077402.3:1..197666,gap(50000),
            NT_187170.1:1..40302,gap(50000),NT_077912.2:1..188020,gap(50000),
            NT_032977.10:1..121390471,gap(50000),NT_187171.1:1..198076,
            gap(100),NT_187172.1:1..278512,gap(100),NT_187173.1:1..2282185,
            gap(100),NT_187174.1:1..63597,gap(100),NT_187175.1:1..83495,
            gap(100),NT_187176.1:1..251763,gap(18000000),
            NT_004487.20:1..80374348,gap(50000),NT_167186.2:1..25337487,
            gap(10000))
ORIGIN"#;
        let contig = contig_text(s.as_bytes()).unwrap();
        assert_eq!(contig.0, &b"ORIGIN"[..]);
        println!("{:?}", contig.1);
    }
}
