use std::str;

use nom::branch::alt;
use nom::bytes::streaming::{is_a, is_not, tag};
use nom::character::complete::none_of;
use nom::character::streaming::{line_ending, one_of, space1};
use nom::combinator::{map, map_res, peek, value};
use nom::error::ErrorKind::MapRes;
use nom::multi::many0;
use nom::sequence::{delimited, pair, preceded, terminated};
use nom::Err::{Error, Incomplete};
use nom::{IResult, Needed, Parser};

use crate::reader::errors::NomParserError;
use crate::reader::location::location;
use crate::seq::{Feature, Location};
use crate::{FeatureKind, QualifierKey};
use crate::reader::field::space_indent;

pub fn features_header(input: &[u8]) -> IResult<&[u8], ()> {
    value(
        (),
        (
            tag("FEATURES"),
            space1,
            tag("Location/Qualifiers"),
            line_ending,
        ),
    ).parse(input)
}

/// Used by `qualifier_value_bare` and `contig`
fn qualifier_value_bare_bytes(indent: usize) -> impl FnMut(&[u8]) -> IResult<&[u8], Vec<u8>> {
    move |mut i| {
        let mut res = Vec::with_capacity(300);
        loop {
            let mut consumed = 0;
            for &b in i {
                match b {
                    b'\r' | b'\n' => break,
                    x => {
                        res.push(x);
                        consumed += 1;
                    }
                }
            }
            i = &i[consumed..];
            if i.is_empty() {
                return Err(Incomplete(Needed::new(1)));
            }
            let (j, _) = line_ending(i)?;
            i = j;
            // check if the next line isn't indented or starts with '/', stop if so
            match preceded(space_indent(indent), peek(none_of("/"))).parse(i) {
                Ok((j, _)) => {
                    i = j;
                }
                Err(Incomplete(n)) => return Err(Incomplete(n)),
                Err(_) => break,
            }
        }
        Ok((i, res))
    }
}

fn location_text<'a>(indent: usize) -> impl Parser<&'a [u8], Output = Location, Error = nom::error::Error<&'a [u8]>> {
    map_res(qualifier_value_bare_bytes(indent), |r| match location(&r) {
        Ok((_, p)) => Ok(p),
        Err(_) => Err(NomParserError::Location),
    })
}

// an unquoted value, may not contain ", but may span multiple lines
// in some cases.
// see: http://www.insdc.org/documents/feature_table.html
fn qualifier_value_bare(indent: usize) -> impl FnMut(&[u8]) -> IResult<&[u8], String> {
    move |i| {
        let (i, mut s) = map_res(qualifier_value_bare_bytes(indent), String::from_utf8).parse(i)?;
        s.shrink_to_fit();
        Ok((i, s))
    }
}

fn qualifier_value_quoted_escape<'a>(
    indent: usize,
) -> impl Parser<&'a [u8], Output = u8, Error = nom::error::Error<&'a [u8]>> {
    alt((
        value(b'"', tag("\"\"")),
        // "" split over a line boundary (does this even happen?)
        value(
            b'"',
            delimited(
                tag("\""),
                pair(line_ending, space_indent(indent)),
                tag("\""),
            ),
        ),
        value(b'\n', pair(line_ending, space_indent(indent))),
    ))
}

// From: http://www.insdc.org/files/feature_table.html#3.3.3
// Since qualifiers convey many different types of information, there are several value formats:
//     1. Free text
//     2. Controlled vocabulary or enumerated values
//     3. Citation or reference numbers
//     4. Sequences

// Just for parsing, to make sure we treat their values correctly if possible
enum QualifierKind {
    FreeText,
    Sequence,
}

impl<'a> From<&'a QualifierKey> for QualifierKind {
    fn from(q: &'a QualifierKey) -> QualifierKind {
        match *q {
            qualifier_key!("translation") => QualifierKind::Sequence,
            _ => QualifierKind::FreeText,
        }
    }
}

// a quoted value, may span multiple lines, quotes
// escaped with ""
fn qualifier_value_quoted<'a, 'b>(
    indent: usize,
    kind: &'a QualifierKind,
) -> impl FnMut(&'b [u8]) -> IResult<&'b [u8], String> + 'a {
    move |input| {
        let (mut i, _) = tag("\"")(input)?;
        let mut res = Vec::with_capacity(2048);
        loop {
            let mut consumed = 0;
            for &b in i {
                match b {
                    b'"' | b'\r' | b'\n' => break,
                    x => {
                        res.push(x);
                        consumed += 1;
                    }
                }
            }
            i = &i[consumed..];
            // try to match an escape sequence
            match qualifier_value_quoted_escape(indent).parse(i) {
                Ok((j, b)) => {
                    // Hacky way to avoid adding linebreaks where they don't belong
                    if b == b'\n' {
                        match *kind {
                            QualifierKind::Sequence => {}
                            _ => res.push(b),
                        }
                    } else {
                        res.push(b);
                    }
                    i = j;
                }
                Err(Incomplete(n)) => return Err(Incomplete(n)),
                Err(_) => break,
            }
        }

        // try to match end of qualifier
        let (j, _) = pair(tag("\""), line_ending).parse(i)?;
        let mut s =
            String::from_utf8(res).map_err(|_| Error(nom::error::Error::new(input, MapRes)))?;
        s.shrink_to_fit(); // save some memory
        Ok((j, s))
    }
}

// needs to be done as a macro because qualifier_key! is a macro
macro_rules! match_qk {
    ($tag:tt) => {
        terminated(value(qualifier_key!($tag), tag($tag)), peek(one_of("=\r\n")))
    };
}

fn qualifier_key(input: &[u8]) -> IResult<&[u8], QualifierKey> {
    alt((
        match_qk!("db_xref"),
        match_qk!("note"),
        match_qk!("gene_synonym"),
        match_qk!("gene"),
        match_qk!("product"),
        match_qk!("transcript_id"),
        match_qk!("codon_start"),
        match_qk!("translation"),
        match_qk!("protein_id"),
        match_qk!("ncRNA_class"),
        match_qk!("pseudo"),
        match_qk!("exception"),
        match_qk!("inference"),
        match_qk!("locus_tag"),
        match_qk!("transl_table"),
        match_qk!("function"),
        map(map_res(is_not("=\r\n"), str::from_utf8), QualifierKey::from),
    )).parse(input)
}

fn qualifier(indent: usize) -> impl FnMut(&[u8]) -> IResult<&[u8], (QualifierKey, Option<String>)> {
    move |i| {
        let (i, _) = space_indent(indent).parse(i)?;
        let (i, _) = tag("/")(i)?;
        let (i, key) = qualifier_key(i)?;
        let (i, val) = alt((
            preceded(
                tag("="),
                map(
                    alt((
                        qualifier_value_quoted(indent, &QualifierKind::from(&key)),
                        qualifier_value_bare(indent),
                    )),
                    Some,
                ),
            ),
            value(None, line_ending),
        )).parse(i)?;
        Ok((i, (key, val)))
    }
}

// needs to be done as a macro because feature_kind! is a macro
macro_rules! match_fk {
    ($tag:tt) => {
        terminated(value(feature_kind!($tag), tag($tag)), peek(one_of(" ")))
    };
}

fn feature_kind(input: &[u8]) -> IResult<&[u8], FeatureKind> {
    alt((
        match_fk!("gene"),
        match_fk!("CDS"),
        match_fk!("mRNA"),
        match_fk!("repeat_region"),
        match_fk!("misc_RNA"),
        match_fk!("ncRNA"),
        map(map_res(is_not(" "), str::from_utf8), FeatureKind::from),
    )).parse(input)
}

pub fn feature(i: &[u8]) -> IResult<&[u8], Feature> {
    let (i, spaces_before) = map(is_a(" "), <[_]>::len).parse(i)?;
    let (i, kind) = feature_kind(i)?;
    let (i, spaces_after) = map(is_a(" "), <[_]>::len).parse(i)?;
    let indent = spaces_before + kind.len() + spaces_after;
    let (i, location) = location_text(indent).parse(i)?;
    let (i, qualifiers) = many0(qualifier(indent)).parse(i)?;
    Ok((
        i,
        Feature {
            kind,
            location,
            qualifiers,
        },
    ))
}

pub fn features(input: &[u8]) -> IResult<&[u8], Vec<Feature>> {
    preceded(features_header, many0(feature)).parse(input)
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_feature_qualifier() {
        let tests = [
            (&b" /test\n"[..], "test", None),
            (&b" /foo=bar\n"[..], "foo", Some("bar")),
            (&b" /foo=\"bar\"\n"[..], "foo", Some("bar")),
            (&b" /foo=\"bar\n bar\"\n"[..], "foo", Some("bar\nbar")),
            (&b" /foo=\"\"\"\"\n"[..], "foo", Some("\"")),
            (&b" /foo=\"test\"\"ing\"\n"[..], "foo", Some("test\"ing")),
            (&b" /foo=\"test\"\n \"ing\"\n"[..], "foo", Some("test\"ing")), // Is this right?
            (
                &b" /rpt_type=long_terminal_repeat\n"[..],
                "rpt_type",
                Some("long_terminal_repeat"),
            ),
        ];
        // Append this, otherwise the parser can't know if the value is complete
        let continuation = b" /continuation";
        for t in &tests {
            println!("Testing: {}", str::from_utf8(t.0).unwrap());
            let v: Vec<u8> = t.0.iter().chain(continuation).cloned().collect();
            match qualifier(1)(&v) {
                Ok((i, q)) => {
                    println!("Got: {:?}", q);
                    assert_eq!(
                        i.len(),
                        continuation.len(),
                        "Complete input not consumed: {}",
                        str::from_utf8(i).unwrap()
                    );
                    assert_eq!(&q.0, t.1);
                    assert_eq!(q.1, t.2.map(String::from));
                }
                Err(x) => {
                    panic!("Parse error: {:?}", x);
                }
            };
        }
    }

    #[test]
    fn test_feature1() {
        test_feature(
            r#"     CDS             join(14924750,14960684..14960875,15034749..15034885,
                     15043989..15044159,15056091..15056280,15060172..15060302,
                     15063572..15063622,15065630..15065797)
                     /gene="KAZN"
                     /gene_synonym="KAZ"
                     /note="Derived by automated computational analysis using
                     gene prediction method: Gnomon."
                     /codon_start=1
                     /product="kazrin isoform X11"
                     /protein_id="XP_016856261.1"
                     /db_xref="GeneID:23254"
                     /db_xref="HGNC:HGNC:29173"
                     /db_xref="HPRD:11122"
                     /translation="MLLREEVSRLQEEVHLLRQMKEMLAKDLEESQGGKSSEVLSATE
                     LRVQLAQKEQELARAKEALQAMKADRKRLKGEKTDLVSQMQQLYATLESREEQLRDFI
                     RNYEQHRKESEDAVKALAKEKDLLEREKWELRRQAKEATDHATALRSQLDLKDNRMKE
                     LEAELAMAKQSLATLTKDVPKRHSLAMPGETVLNGNQEWVVQADLPLTAAIRQSQQTL
                     YHSHPPHPADRQAVRVSPCHSRQPSVISDASAAEGDRSSTPSDINSPRHRTHSLCNGD
                     SPGPVQKNLHNPIVQSLEDLEDQKRKKKKEKMGFGSISRVFARGKQRKSLDPGLFDGT
                     APDYYIEEDADW""#,
        );
    }

    #[test]
    fn test_feature2() {
        test_feature(
            r#"     prim_transcript 44..>579
                     /gene="kin2""#,
        );
    }

    #[test]
    fn test_feature3() {
        test_feature(
            r#"     gene            <1..>209
                     /gene="TAS1R2"
                     /pseudogene="unprocessed""#,
        );
    }

    fn test_feature(f: &str) {
        let f = format!("{}\n     CDS             ", f); // Needs the start of the next feature to know it's finished
        let f = f.as_bytes();
        let f_parsed = feature(f);
        match f_parsed {
            Ok((i, o)) => {
                println!("[{}] => {:?}", str::from_utf8(i).unwrap(), o);
                assert_eq!(i, b"     CDS             ");
            }
            Err(Incomplete(_)) => {
                panic!();
            }
            Err(e) => {
                panic!("{:?}", e);
            }
        }
        incomplete_test(f, b"     CDS             ", &feature);
    }

    fn incomplete_test<T: ::std::fmt::Debug>(
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
                Err(Incomplete(_)) => {}
                x => panic!("{:?} => {:#?}", &input[..n], x),
            }
        }
    }
}
