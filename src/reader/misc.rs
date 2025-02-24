use crate::reader::field::{contig_text, field_string, fields, fill_seq_fields};
use nom::branch::alt;
use nom::bytes::streaming::tag;
use nom::character::streaming::{alpha1, digit1, line_ending, multispace0, multispace1, not_line_ending, space0, space1};
use nom::combinator::{complete, map, map_res, not, opt, value};
use nom::multi::{fold_many0, many0, many1};
use nom::sequence::{delimited, pair, preceded, separated_pair, tuple};
use nom::{IResult, Parser};
use nom::error::ParseError;
use std::{cmp, str};
use crate::reader::feature_table::features;
use crate::reader::locus::locus;
use crate::seq::{REASONABLE_SEQ_LEN, Seq};

pub fn ignored_line(input: &[u8]) -> IResult<&[u8], &str> {
    delimited(
        not(alt((tag("ORIGIN"), tag("CONTIG"), tag("FEATURES")))),
        map_res(not_line_ending, str::from_utf8),
        line_ending,
    ).parse(input)
}

// Added this to skip header info, for example in genbank .seq files. Not sure
// if it's the best/most correct solution
pub fn skip_preamble(i: &[u8]) -> IResult<&[u8], ()> {
    fold_many0(
        delimited(
            not(tag("LOCUS")),
            map_res(not_line_ending, str::from_utf8),
            line_ending,
        ),
        || (),
        |_, i| warn!("Ignoring line: {}", i),
    ).parse(i)
}

pub fn origin_tag(i: &[u8]) -> IResult<&[u8], Option<String>> {
    alt((
        value(None, tuple((tag("ORIGIN"), space0, line_ending))),
        map(field_string(0, "ORIGIN", true), Some),
    )).parse(i)
}

pub fn double_slash(i: &[u8]) -> IResult<&[u8], ()> {
    value((), tag("//")).parse(i)
}

// ORIGIN

fn origin<'a>(len: Option<usize>) -> impl Parser<&'a [u8], Output = (Option<String>, Vec<u8>), Error = nom::error::Error<&'a [u8]>> {
    separated_pair(
        origin_tag,
        space1,
        sequence(len)
    )
}

fn sequence_chunk(i: &[u8]) -> IResult<&[u8], &[u8]> {
    delimited(opt(pair(digit1, space1)), alpha1, multispace1).parse(i)
}

fn sequence<'a>(len: Option<usize>) -> impl Parser<&'a [u8], Output = Vec<u8>, Error = nom::error::Error<&'a [u8]>> {
    map_res(
        fold_many0(sequence_chunk,
                   move || match len {
                       Some(len) => {
                           Vec::with_capacity(cmp::min(len, REASONABLE_SEQ_LEN))
                       },
                       None => Vec::new()
                   },
                   |mut acc: Vec<u8>, item: &[u8]| {
                       acc.extend_from_slice(item);
                       acc
                   }
        ),
        move |res: Vec<u8>| -> Result<Vec<u8>, ()> {
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
}

// I think this is deprecated but sometimes it comes after features,
// let's just ignore it.
pub fn base_count(i: &[u8]) -> IResult<&[u8], ()> {
    value((), field_string(0, "BASE COUNT", false)).parse(i)
}

fn gb(i: &[u8]) -> IResult<&[u8], Seq> {
    let (i, locus) =  locus(i)?;
    let (i, metadata) =  map_res(fields, |f| fill_seq_fields(Seq::empty(), f)).parse(i)?;
    let (i, _) = many0(ignored_line).parse(i)?;
    let (i, features) =  opt(features).parse(i)?;
    let (i, _) = opt(base_count).parse(i)?;
    let (i, contig) =  opt(contig_text).parse(i)?;
    let (i, origin) =  opt(origin(locus.len)).parse(i)?;
    let (i, _) = multispace0(i)?;
    let (i, _) = tag("//")(i)?;
    let (i, _) = opt(complete(multispace1)).parse(i)?;
    Ok((i, Seq {
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
    }))
}

pub fn gb_records(i: &[u8]) -> IResult<&[u8], Vec<Seq>> {
    preceded(skip_preamble, many1(gb)).parse(i)
}