use itertools::Itertools;
use nom;
use nom::types::CompleteByteSlice;
use nom::{alpha, digit, line_ending, multispace, not_line_ending, space, IResult};
use crate::seq::{
    After, Before, Date, Feature, FeatureKind, Location, QualifierKey, Reference, Seq, Source, Topology,
    REASONABLE_SEQ_LEN, GapLength
};
use std::cmp;
use std::iter::once;
use std::str;

use crate::reader::errors::NomParserError;

// A few helper functions, used throughout

named!(
    numeric_usize<usize>,
    flat_map!(recognize!(digit), parse_to!(usize))
);

macro_rules! numeric_i64 (
    ($i:expr, $(,)*) => (
        flat_map!(
          $i,
          recognize!(do_parse!(opt!(tag!("-")) >> digit >> ())),
          parse_to!(i64)
        )
    )
);

/// Matches `indent` spaces
fn space_indent(i: &[u8], indent: usize) -> IResult<&[u8], ()> {
    if indent == 0 {
        return Ok((i, ()));
    }
    let mut consumed = 0;
    for &b in i {
        match b {
            b' ' => {
                consumed += 1;
                if consumed == indent {
                    return Ok((&i[consumed..], ()));
                }
            }
            _ => return Err(nom::Err::Error(nom::Context::Code(i, nom::ErrorKind::Tag))),
        }
    }
    Err(nom::Err::Incomplete(nom::Needed::Unknown))
}

// LOCUS

/// Internal type for convenience, we move everything into the `Seq` struct later
#[derive(Debug)]
pub struct Locus {
    pub name: Option<String>,
    pub len: Option<usize>,
    pub topology: Topology,
    pub date: Option<Date>,
    pub molecule_type: Option<String>,
    pub division: String,
}

named!(
    topology<Topology>,
    alt!(
        tag!("linear")   => { |_| Topology::Linear   } |
        tag!("circular") => { |_| Topology::Circular }
)
);

named!(
    gb_month<usize>,
    alt!(
        tag!("JAN") => { |_| 1  } |
        tag!("FEB") => { |_| 2  } |
        tag!("MAR") => { |_| 3  } |
        tag!("APR") => { |_| 4  } |
        tag!("MAY") => { |_| 5  } |
        tag!("JUN") => { |_| 6  } |
        tag!("JUL") => { |_| 7  } |
        tag!("AUG") => { |_| 8  } |
        tag!("SEP") => { |_| 9  } |
        tag!("OCT") => { |_| 10 } |
        tag!("NOV") => { |_| 11 } |
        tag!("DEC") => { |_| 12 }
    )
);

named!(
    date<Date>,
    add_return_error!(
        ErrorKind::Custom(NomParserError::Date as u32),
        map_res!(
            do_parse!(
                day: numeric_usize
                    >> tag!("-")
                    >> month: gb_month
                    >> tag!("-")
                    >> year: numeric_i64!()
                    >> (Date::from_ymd(year as i32, month as u32, day as u32))
            ),
            |res| res
        )
    )
);

named!(locus_name<&str>, to_str!(is_not!(": \t\r\n")));

named!(spaces, eat_separator!(&b" \t"[..]));

named!(molecule_type<&str>, to_str!(is_not!(" ")));

named!(
    locus_full<Locus>,
    sep!(
        spaces,
        do_parse!(
            name: locus_name
                >> len: numeric_usize
                >> tag!("bp")
                >> mol: molecule_type
                >> topology: topology
                >> section: opt!(to_str!(alpha))
                >> date: opt!(date)
                >> (Locus {
                    name: Some(String::from(name)),
                    len: Some(len),
                    topology,
                    date,
                    molecule_type: Some(mol.into()),
                    division: section.unwrap_or("UNK").into(),
                })
        )
    )
);

// EMBL style LOCUS, no topology info
named!(
    locus_traditional<Locus>,
    sep!(
        spaces,
        do_parse!(
            name: locus_name >>
            len: numeric_usize >>
            tag!("bp") >>
            mol: molecule_type >> // TODO: what are other molecule types?
            section: to_str!(alpha) >> date: date >> (Locus {
                name: Some(String::from(name)),
                len: Some(len),
                topology: Topology::Linear,
                date: Some(date),
                molecule_type: Some(mol.into()),
                division: section.into(),
            })
        )
    )
);

// Just give up :)
named!(
    locus_tag_only<Locus>,
    do_parse!(
        stuff: opt!(to_str!(not_line_ending))
            >> ({
                warn!("Failed to parse: {:?}", stuff.unwrap_or_default()); //TODO: fix
                Locus {
                    name: None,
                    len: None,
                    topology: Topology::Linear,
                    date: None,
                    molecule_type: None,
                    division: "UNK".into(),
                }
            })
    )
);

// Added this to skip header info, for example in genbank .seq files. Not sure
// if it's the best/most correct solution
named!(pub skip_preamble<()>,
    fold_many0!(
        do_parse!(
            not!(tag!("LOCUS")) >>
                ignored_line: to_str!(not_line_ending) >>
                line_ending >>
                (ignored_line)
        ), (),
        | _, i | warn!("Ignoring line: {}", i)
    )
);

named!(
    pub locus<Locus>,
    do_parse!(
        tag!("LOCUS") >> space
        >> metadata: alt!(locus_full | locus_traditional | locus_tag_only)
        >> line_ending >> (metadata)
    )
);

// Metadata

named_args!(field<'a>(indent: usize, name: &str, keep_ws: bool) <String>,
    to_string!(apply!(field_bytes, indent, name, keep_ws))
);

/// Concatenates lines, optionally interpolating with newlines
fn concat_lines<'a, T: Iterator<Item = &'a [u8]>>(lines: T, keep_linebreaks: bool) -> Vec<u8> {
    if keep_linebreaks {
        Itertools::intersperse(lines, b"\n")
            .flatten()
            .cloned()
            .collect()
    } else {
        lines.flatten().cloned().collect()
    }
}

// Matches an entry beginning with an (optionally) indented identifier followed
// by one of more indented lines of free text
named_args!(field_bytes<'a>(indent: usize, name: &str, keep_ws: bool) <Vec<u8>>,
    do_parse!(
        apply!(space_indent, indent) >>
        tag!(name) >>
        spaces: map!(is_a!(" "), <[_]>::len) >>
        first_line: not_line_ending >>
        line_ending >>
        more_lines: many0!(
            delimited!(
                apply!(space_indent, indent + spaces + name.len()),
                not_line_ending,
                line_ending
                )
            ) >>
        (concat_lines(once(first_line).chain(more_lines.into_iter()), keep_ws))
    )
);

// Just for convenience
named_args!(toplevel_field<'a>(name: &str, keep_ws: bool) <String>,
    apply!(field, 0, name, keep_ws)
);

// REFERENCE

named!(
    reference<Reference>,
    do_parse!(
        description: apply!(field, 0, "REFERENCE", true)
            >> authors: opt!(apply!(field, 2, "AUTHORS", true))
            >> consortium: opt!(apply!(field, 2, "CONSRTM", true))
            >> title: apply!(field, 2, "TITLE", true)
            >> journal: opt!(apply!(field, 2, "JOURNAL", true))
            >> pubmed: opt!(apply!(field, 3, "PUBMED", false))
            >> remark: opt!(apply!(field, 2, "REMARK", true))
            >> (Reference {
                description,
                authors,
                consortium,
                title,
                journal,
                pubmed,
                remark,
            })
    )
);

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
    ($i:expr, $field: ident, $keep_ws:expr) => {
        do_parse!(
            $i,
            s: apply!(toplevel_field, stringify!($field), $keep_ws) >> (Field::$field(s))
        )
    };
}

named!(
    source<Field>,
    do_parse!(
        source: apply!(field, 0, "SOURCE", true)
            >> organism: opt!(apply!(field, 2, "ORGANISM", true))
            >> (Field::SOURCE(Source {
                source,
                organism,
            }))
    )
);

named!(
    pub any_field<Field>,
    alt!(
        parse_field!(DEFINITION, true) |
        parse_field!(ACCESSION, true)  |
        parse_field!(VERSION, true)    |
        parse_field!(DBLINK, true)     |
        parse_field!(KEYWORDS, true)   |
        source |
        call!(reference) => { Field::REFERENCE } |
        parse_field!(COMMENT, true) |
        //TODO: unrecognised lines?
        do_parse!(line: ignored_line >>
                    (Field::UnrecognisedLine(line.into())))
    )
);

named!(fields<Vec<Field>>, many0!(any_field));

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

// Features table

named!(
    pub features_header<()>,
    do_parse!(tag!("FEATURES") >> multispace >> tag!("Location/Qualifiers") >> line_ending >> ())
);

named!(
    ignored_line<&str>,
    do_parse!(
        not!(alt!(tag!("ORIGIN") | tag!("CONTIG") | tag!("FEATURES"))) >> // TODO: fix this
        content: to_str!(not_line_ending) >> line_ending >> (content)
    )
);

/// Used by `qualifier_value_bare` and `contig`
fn qualifier_value_bare_bytes(input: &[u8], indent: usize) -> IResult<&[u8], Vec<u8>> {
    let mut i = input;
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
            return Err(nom::Err::Incomplete(nom::Needed::Size(1)));
        }
        let (j, _) = line_ending(i)?;
        i = j;
        // check if the next line isn't indented or starts with '/', stop if so
        match do_parse!(
            i,
            apply!(space_indent, indent) >> peek!(none_of!("/")) >> ()
        ) {
            Ok((j, _)) => {
                i = j;
            }
            Err(nom::Err::Incomplete(n)) => return Err(nom::Err::Incomplete(n)),
            Err(_) => break,
        }
    }
    Ok((i, res))
}

fn pos_text(i: &[u8], indent: usize) -> IResult<&[u8], Location> {
    map_res!(
        i,
        apply!(qualifier_value_bare_bytes, indent),
        |r: Vec<u8>| match location(CompleteByteSlice(&r)) {
            Ok((_, p)) => Ok(p),               // TODO: Check we're matching the whole input
            Err(e) => Err(format!("{:?}", e)), // TODO: Pass the error on somehow
        }
    )
}

// an unquoted value, may not contain ", but may span multiple lines
// in some cases.
// see: http://www.insdc.org/documents/feature_table.html
fn qualifier_value_bare(i: &[u8], indent: usize) -> IResult<&[u8], String> {
    let (i, mut s) = to_string!(i, apply!(qualifier_value_bare_bytes, indent))?;
    s.shrink_to_fit();
    Ok((i, s))
}

named_args!(qualifier_value_quoted_escape<'a>(indent: usize) <u8>,
    alt!(
        tag!("\"\"") => { |_| b'"' } |
        // "" split over a line boundary (does this even happen?)
        delimited!(
            tag!("\""),
            pair!(line_ending, apply!(space_indent, indent)),
            tag!("\"")
        ) => { |_| b'"' } |
        pair!(line_ending, apply!(space_indent, indent)) => {|_| b'\n' }
    )
);

// a quoted value, may span multiple lines, quotes
// escaped with ""
fn qualifier_value_quoted<'a>(
    input: &'a [u8],
    indent: usize,
    kind: &QualifierKind,
) -> IResult<&'a [u8], String> {
    let (mut i, _) = tag!(input, "\"")?;
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
        match qualifier_value_quoted_escape(i, indent) {
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
            Err(nom::Err::Incomplete(n)) => return Err(nom::Err::Incomplete(n)),
            Err(_) => break,
        }
    }

    // try to match end of qualifier
    let (j, _) = pair!(i, tag!("\""), line_ending)?;
    // Convert to String like this to get the right error (with context)
    let (_, mut s) = to_string!(input, value!(res))?;
    s.shrink_to_fit(); // save some memory
    Ok((j, s))
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

macro_rules! match_tags (
    ($i:expr, $atom:ident!, $term:ident!($($term_args:tt)*), $first:tt | $($rest:tt)*) => (
        alt!($i, terminated!(tag!($first), peek!($term!($($term_args)*))) =>
             { |_| $atom!($first) } | match_tags!($atom!, $term!($($term_args)*), $($rest)*))
    );
    ($i:expr, $atom:ident!, $term:ident!($($term_args:tt)*), $macro:ident!($($args:tt)*)) => (
        alt!($i, $macro!($($args)*))
    )
);

named!(
    qualifier_key<QualifierKey>,
    match_tags!(qualifier_key!,
      one_of!("=\r\n"),
      "db_xref"       |
      "note"          |
      "gene_synonym"  |
      "gene"          |
      "product"       |
      "transcript_id" |
      "codon_start"   |
      "translation"   |
      "protein_id"    |
      "ncRNA_class"   |
      "pseudo"        |
      "exception"     |
      "inference"     |
      "locus_tag"     |
      "transl_table"  |
      "function"      |
      map!(to_str!(is_not!("=\r\n")), QualifierKey::from)
    )
);

fn qualifier(i: &[u8], indent: usize) -> IResult<&[u8], (QualifierKey, Option<String>)> {
    do_parse!(
        i,
        apply!(space_indent, indent)
            >> tag!("/")
            >> key: qualifier_key
            >> val: alt!(
                    preceded!(
                        tag!("="),
                        alt!(
                            apply!(qualifier_value_quoted, indent, &QualifierKind::from(&key))  |
                            apply!(qualifier_value_bare, indent)
                            )
                    ) => { Some } |
                    line_ending => {|_| None }
              )
            >> (key, val)
    )
}

fn qualifiers(
    i: &[u8],
    indent: usize,
) -> IResult<&[u8], Vec<(QualifierKey, Option<String>)>> {
    fold_many0!(
        i,
        apply!(qualifier, indent),
        Vec::new(),
        |mut acc: Vec<(QualifierKey, Option<String>)>, item: (QualifierKey, Option<String>)| {
            acc.push(item);
            acc
        }
    )
}

named!(
    feature_kind<FeatureKind>,
    match_tags!(
           feature_kind!,
           tag!(" "),
           "gene" | "CDS" | "mRNA" | "repeat_region" | "misc_RNA" | "ncRNA" |
           map!(to_str!(is_not!(" ")), FeatureKind::from)
      )
);

named!(
    pub feature <Feature>,
    do_parse!(
        spaces_before: map!(is_a!(" "), <[_]>::len) >> kind: call!(feature_kind)
            >> spaces_after: map!(is_a!(" "), <[_]>::len)
            >> indent: value!(spaces_before + kind.len() + spaces_after)
            >> location: apply!(pos_text, indent)
            >> qualifiers: apply!(qualifiers, indent) >> (Feature {
            kind,
            location,
            qualifiers,
        })
    )
);

named!(
    features<Vec<Feature>>,
    do_parse!(features_header >> features: many0!(feature) >> (features))
);

// Feature locations

named!(
    pos_single<CompleteByteSlice, Location>,
    map!(numeric_i64!(), |i| Location::Range(((i - 1), Before(false)), (i, After(false)))) // Convert from 1-based
); // to 0-based

named!(
    pos_between<CompleteByteSlice, Location>,
    map_res!(
        separated_pair!(numeric_i64!(), tag!("^"), numeric_i64!()),
        |(a, b): (i64, i64)| -> Result<Location, String> {
            // This check is not 100% foolproof, because if b == 1, a must be
            // the last nt of a circular sequence, but since we might not yet
            // know the sequence length, we can't check this yet
            if (a - b).abs() == 1 || ((a == 1) ^ (b == 1))
            {
                Ok(Location::Between(a - 1, b - 1))
            } else {
                Err(format!(
                    "Invalid location, coordinates separated by ^ must be adjacent"
                ))
            }
        }
    )
);

// A range x..y
named!(
    pos_span<CompleteByteSlice, Location>,
    do_parse!(
        before: opt!(char!('<')) >>
        a: numeric_i64!() >>
        tag!("..") >>
        after: opt!(char!('>')) >>
        b: numeric_i64!() >>
        (Location::Range((a - 1, Before(before.is_some())), (b, After(after.is_some()))))
    )
);

macro_rules! compound_pos {
    ($i:expr, $kind: ident, $tag:expr) => {
        do_parse!(
            $i,
            tag!($tag)
                >> locations:
                    delimited!(tag!("("), separated_list!(tag!(","), location), tag!(")"))
                >> (Location::$kind(locations))
        )
    };
}

named!(pos_join<CompleteByteSlice, Location>, compound_pos!(Join, "join"));
named!(pos_order<CompleteByteSlice, Location>, compound_pos!(Order, "order"));
named!(pos_bond<CompleteByteSlice, Location>, compound_pos!(Bond, "bond"));
named!(pos_oneof<CompleteByteSlice, Location>, compound_pos!(OneOf, "one-of"));

named!(
    pos_complement<CompleteByteSlice, Location>,
    do_parse!(
        tag!("complement") >> location: delimited!(tag!("("), location, tag!(")"))
            >> (Location::Complement(Box::new(location)))
    )
);

named!(
    pos_external<CompleteByteSlice, Location>,
    do_parse!(
        accession: to_str!(map!(is_not!(": \t\r\n"), |x| x.0)) // convert to &[u8]
            >> location:
                opt!(preceded!(
                    tag!(":"),
                    alt_complete!(
                        // everything except pos_external
                        pos_span | pos_join | pos_complement |
                        pos_between | pos_single | pos_gap | 
                        pos_order | pos_oneof | pos_bond
                    )
                )) >> (Location::External(accession.to_owned(), location.map(Box::new)))
    )
);

// This is probably only legal in a CONTIG field, but then again, who knows, it
// can't hurt to accept it everywhere
named!(
    pos_gap<CompleteByteSlice, Location>,
    do_parse!(tag!("gap(") >> length: gap_length >> tag!(")") >> (Location::Gap(length)))
);

named!(
    gap_length<CompleteByteSlice, GapLength>,
    alt_complete!(
        numeric_i64!() => { GapLength::Known } |
        tag!("unk100") => { |_| GapLength::Unk100 } |
        tag!("") => {|_| GapLength::Unknown }
    )
);

named!(
    pub location<CompleteByteSlice, Location>,
    // Order is important here, `pos_single` eats the first number in a
    // range otherwise
    alt_complete!(
        pos_span | pos_join | pos_complement | pos_between |
        pos_single |  pos_gap | pos_order | pos_oneof | 
        pos_bond | pos_external
    )
);

// CONTIG

named!(
    pub contig_text<Location>,
    map_res_custom_error!(
        NomParserError::Location,
        apply!(field_bytes, 0, "CONTIG", false),
        |r: Vec<u8>| match location(CompleteByteSlice(&r)) {
            Ok((_, p)) => Ok(p), // TODO: Check we're matching the whole input
            Err(e) => Err(format!("{:?}", e)) // TODO: Pass this error on once
                                              // we have proper custom errors
        }
    )
);

// ORIGIN

named!(
    pub origin_tag<Option<String> >,
    alt!(
        do_parse!(
            tag!("ORIGIN") >>
            opt!(space) >>
            line_ending >>
            (None)
        ) |
        apply!(toplevel_field, "ORIGIN", true) => { Some }
    )
);

named_args!(origin(len: Option<usize>)<(Option<String>, Vec<u8>)>,
       do_parse!(
               origin_text: origin_tag >>
               space >>
               seq: apply!(sequence, len) >>
               (origin_text, seq)
       )
);

named!(
    sequence_chunk,
    delimited!(opt!(pair!(digit, space)), alpha, multispace)
);

named_args!(sequence(len: Option<usize>) <Vec<u8>>,
    map_res_custom_error!(
        NomParserError::SequenceLength,
        fold_many0!(sequence_chunk,
            match len {
                Some(len) => {
                    Vec::with_capacity(cmp::min(len, REASONABLE_SEQ_LEN))
                },
                None => Vec::new()
            },
            | mut acc: Vec<u8>, item: &[u8] | {
                acc.extend_from_slice(item);
                acc
            }
        ),
        |res: Vec<u8>| -> Result<Vec<u8>, ()> {
            if let Some(len) = len {
                if res.len() != len {
                    Err(())
                } else {
                    Ok(res)
                }
            } else {
                Ok(res)
            }
        }
    )
);

// I think this is deprecated but sometimes it comes after features,
// let's just ignore it.
named!(pub base_count<()>,
       value!((), apply!(toplevel_field, "BASE COUNT", false))
);

named!(
    gb<Seq>,
    do_parse!(
        locus: locus
            >> metadata: map_res!(fields, |f| fill_seq_fields(Seq::empty(), f))
            >> many0!(ignored_line)
            >> features: opt!(features)
            >> opt!(call!(base_count))
            >> contig: opt!(contig_text)
            >> origin: opt!(apply!(origin, locus.len))
            >> opt!(multispace)
            >> tag!("//")
            >> opt!(complete!(multispace))
            >> (Seq {
                name: locus.name,
                date: locus.date,
                topology: locus.topology,
                len: locus.len,
                molecule_type: locus.molecule_type,
                division: locus.division,
                seq: origin.map(|o| o.1).unwrap_or_default(),
                contig,
                features: features.unwrap_or_default(),
                ..metadata
            })
    )
);

named!(pub gb_records<Vec<Seq>>, preceded!(skip_preamble, many1!(gb)));

named!(pub line_ending_type_hack<()>, value!((), line_ending));
named!(pub double_slash<()>, value!((), tag!("//")));

#[cfg(test)]
mod test {
    use super::*;
    #[test]
    fn test_date() {
        let d = date(&b"01-AUG-2014\n"[..]);
        assert_eq!(d, Ok((&b"\n"[..], Date::from_ymd(2014, 8, 1).unwrap())))
    }

    #[test]
    fn test_gap() {
        assert_eq!(location(CompleteByteSlice(b"gap()")), 
            Ok((CompleteByteSlice(b""), Location::Gap(GapLength::Unknown))));
        assert_eq!(location(CompleteByteSlice(b"gap(10)")), 
            Ok((CompleteByteSlice(b""), Location::Gap(GapLength::Known(10)))));
        assert_eq!(location(CompleteByteSlice(b"gap(unk100)")), 
            Ok((CompleteByteSlice(b""), Location::Gap(GapLength::Unk100))));
    }

    #[test]
    fn test_parse_locus() {
        let loci = [
            "LOCUS       U00096               4641652 bp    DNA     circular BCT 01-AUG-2014\n",
            "LOCUS       U00096               4641652 bp    DNA     circular 01-AUG-2014\n",
            "LOCUS       pBAD30       4923 bp    DNA   circular            \n",
        ];
        for l in &loci {
            println!("Testing: {}", l);
            let lb = l.as_bytes();
            let loc = locus(lb);
            println!("{:?}", loc);
            match loc {
                Ok((i, _)) => assert!(i.is_empty()),
                Err(e) => {
                    panic!("{:?}", e);
                }
            }
            // incomplete locus test
            incomplete_test(lb, b"", &locus);
        }
        assert_eq!(
            locus(loci[2].as_bytes()).unwrap().1.topology,
            Topology::Circular
        );
    }

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
            match qualifier(&v, 1) {
                Ok((i, q)) => {
                    println!("Got: {:?}", q);
                    assert!(
                        i.len() == continuation.len(),
                        "Complete input not consumed: {}",
                        str::from_utf8(i).unwrap()
                    );
                    assert!(q.0 == QualifierKey::from(t.1));
                    assert!(q.1 == t.2.map(String::from));
                }
                Err(x) => {
                    panic!("Parse error: {:?}", x);
                }
            };
        }
    }

    #[test]
    fn test_feature() {
        let f = r#"     CDS             join(14924750,14960684..14960875,15034749..15034885,
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
                     APDYYIEEDADW"
     CDS             "#; // Needs the start of the next feature to know it's finished
        /* let f = r#"     repeat_region   complement(5293260..5294954)
                     /note="ERV-9; Derived by automated computational analysis
                     using gene prediction method: Curated Genomic."
                     /rpt_type=long_terminal_repeat
     gene            complement(5301014..5301946)
                     /gene="OR51B4"
                     /gene_synonym="HOR5'Beta1"
                     /note="olfactory receptor family 51 subfamily B member 4;
                     Derived by automated computational analysis using gene
                     prediction method: BestRefSeq."
                     /db_xref="GeneID:79339"
                     /db_xref="HGNC:HGNC:14708"
     mRNA            "#; */

        let f = f.as_bytes();
        let f_parsed = feature(f);
        match f_parsed {
            Ok((i, o)) => {
                println!("[{}] => {:?}", str::from_utf8(i).unwrap(), o);
                assert!(i == b"     CDS             ");
            }
            Err(nom::Err::Incomplete(_)) => {
                panic!();
            }
            Err(e) => {
                panic!("{:?}", e);
            }
        }
        incomplete_test(f, b"     CDS             ", &feature);
    }

    #[test]
    fn test_pos_single() {
        let s = "1\n";
        let p = location(CompleteByteSlice(s.as_bytes()));
        match p {
            Ok((_, p)) => match p {
                Location::Range((0, Before(false)), (1, After(false))) => {}
                x => {
                    panic!("{:?}", x);
                }
            },
            x => {
                panic!("{:?}", x);
            }
        }
    }

    #[test]
    fn test_pos_join() {
        let s = "join(1,3)";
        let p = pos_join(CompleteByteSlice(s.as_bytes()));
        let inner = vec![Location::simple_range(0, 1), Location::simple_range(2, 3)];
        assert_eq!(
            p,
            Ok(((CompleteByteSlice(&b""[..])), Location::Join(inner)))
        );
    }

    fn incomplete_test<T: ::std::fmt::Debug>(
        input: &[u8],
        remainder: &[u8],
        parser: impl Fn(&[u8]) -> IResult<&[u8], T>,
    ) {
        match parser(input) {
            Ok((i, _)) => {
                assert!(i == remainder);
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
        incomplete_test(f, b"FEATURES             Location/Qualifiers", &any_field);
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
