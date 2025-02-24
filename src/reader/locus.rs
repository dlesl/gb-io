use std::str;

use nom::branch::alt;
use nom::bytes::streaming::{is_not, tag};
use nom::character::streaming::{
    alpha1, i32, line_ending, not_line_ending, space0, space1, u32, u64,
};
use nom::combinator::{map, map_res, opt, value};
use nom::sequence::{delimited, preceded, tuple};
use nom::{IResult, Parser};

use crate::seq::{Date, Topology};

#[derive(Debug)]
pub struct Locus {
    pub name: Option<String>,
    pub len: Option<usize>,
    pub topology: Topology,
    pub date: Option<Date>,
    pub molecule_type: Option<String>,
    pub division: String,
}

fn topology(input: &[u8]) -> IResult<&[u8], Topology> {
    alt((
        map(tag("linear"), |_| Topology::Linear),
        map(tag("circular"), |_| Topology::Circular),
    )).parse(input)
}

fn date(input: &[u8]) -> IResult<&[u8], Date> {
    let month_parser = alt((
        value(1, tag("JAN")),
        value(2, tag("FEB")),
        value(3, tag("MAR")),
        value(4, tag("APR")),
        value(5, tag("MAY")),
        value(6, tag("JUN")),
        value(7, tag("JUL")),
        value(8, tag("AUG")),
        value(9, tag("SEP")),
        value(10, tag("OCT")),
        value(11, tag("NOV")),
        value(12, tag("DEC")),
    ));

    map_res(
        tuple((u32, tag("-"), month_parser, tag("-"), i32)),
        |(day, _, month, _, year)| Date::from_ymd(year, month, day),
    ).parse(input)
}

fn name(input: &[u8]) -> IResult<&[u8], &str> {
    map_res(is_not(": \t\r\n"), str::from_utf8).parse(input)
}

fn molecule_type(input: &[u8]) -> IResult<&[u8], &str> {
    map_res(is_not(" "), str::from_utf8).parse(input)
}

fn locus_full(input: &[u8]) -> IResult<&[u8], Locus> {
    map(
        tuple((
            name,
            preceded(space1, u64),
            preceded(space1, tag("bp")),
            preceded(space1, molecule_type),
            preceded(space1, topology),
            opt(preceded(space1, map_res(alpha1, str::from_utf8))),
            opt(preceded(space1, date)),
        )),
        |(name, len, _, molecule_type, topology, division, date)| Locus {
            name: Some(String::from(name)),
            len: Some(len as usize),
            topology,
            date,
            molecule_type: Some(molecule_type.into()),
            division: division.unwrap_or("UNK").into(),
        },
    ).parse(input)
}

// EMBL style LOCUS, no topology info
fn locus_traditional(input: &[u8]) -> IResult<&[u8], Locus> {
    map(
        tuple((
            name,
            preceded(space1, u64),
            preceded(space1, tag("bp")),
            preceded(space1, molecule_type),
            map_res(preceded(space1, alpha1), str::from_utf8),
            opt(preceded(space1, date)),
        )),
        |(name, len, _, molecule_type, division, date)| Locus {
            name: Some(String::from(name)),
            len: Some(len as usize),
            topology: Topology::Linear,
            date,
            molecule_type: Some(molecule_type.into()),
            division: division.into(),
        },
    ).parse(input)
}

// Just give up :)
fn locus_tag_only(input: &[u8]) -> IResult<&[u8], Locus> {
    map(
        opt(map(not_line_ending, String::from_utf8_lossy)),
        |stuff| {
            warn!("Failed to parse: {:?}", stuff.unwrap_or_default()); //TODO: fix
            Locus {
                name: None,
                len: None,
                topology: Topology::Linear,
                date: None,
                molecule_type: None,
                division: "UNK".into(),
            }
        },
    ).parse(input)
}

pub fn locus(input: &[u8]) -> IResult<&[u8], Locus> {
    delimited(
        tuple((tag("LOCUS"), space1)),
        alt((locus_full, locus_traditional, locus_tag_only)),
        tuple((space0, line_ending)),
    ).parse(input)
}

#[cfg(test)]
mod test {
    use nom::Err::Incomplete;

    use super::*;

    #[test]
    fn test_date() {
        let d = date(&b"01-AUG-2014\n"[..]);
        assert_eq!(d, Ok((&b"\n"[..], Date::from_ymd(2014, 8, 1).unwrap())))
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
