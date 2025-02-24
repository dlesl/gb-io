/// Parsers for location specifiers
use std::str;

use branch::alt;
use bytes::complete::{is_not, tag};
use character::complete::{char, i64};
use combinator::{map, map_opt, map_res, opt};
use multi::separated_list1;
use nom::combinator::value;
use nom::{IResult, Parser};
use nom::{branch, bytes, character, combinator, multi, sequence};
use sequence::{delimited, preceded, separated_pair, tuple};

use crate::seq::{After, Before, GapLength, Location};

fn location_single(input: &[u8]) -> IResult<&[u8], Location> {
    map(i64, |i| {
        Location::Range(((i - 1), Before(false)), (i, After(false)))
    }).parse(input)
}

fn location_between(input: &[u8]) -> IResult<&[u8], Location> {
    map_opt(
        separated_pair(i64, char('^'), i64),
        // This check is not 100% foolproof, because if b == 1, a must be
        // the last nt of a circular sequence, but since we might not yet
        // know the sequence length, we can't check this yet
        |(a, b)| {
            ((a - b).abs() == 1 || ((a == 1) ^ (b == 1))).then_some(Location::Between(a - 1, b - 1))
        },
    ).parse(input)
}

fn location_span(input: &[u8]) -> IResult<&[u8], Location> {
    let parser = tuple((opt(char('<')), i64, tag(".."), opt(char('>')), i64));
    map(parser, |(before, a, _, after, b)| {
        Location::Range(
            (a - 1, Before(before.is_some())),
            (b, After(after.is_some())),
        )
    }).parse(input)
}

// fn operator<'a, F, O>(name: &'static str, inner: F) -> impl FnMut(&'a [u8]) -> IResult<&'a [u8], O>
// where
//     F: FnMut(&'a [u8]) -> IResult<&'a [u8], O>,
// {
//     preceded(tag(name), delimited(char('('), inner, char(')')))
// }

fn location_compound(input: &[u8]) -> IResult<&[u8], Location> {
    todo!("location_compound")
    // fn parser<'a>(
    //     name: &'static str,
    //     constructor: impl Fn(Vec<Location>) -> Location,
    // ) -> impl FnMut(&'a [u8]) -> IResult<&[u8], Location> {
    //     map(
    //         operator(name, separated_list1(char(','), location)),
    //         constructor,
    //     )
    // }
    // alt((
    //     parser("join", Location::Join),
    //     parser("order", Location::Order),
    //     parser("bond", Location::Bond),
    //     parser("one-of", Location::OneOf),
    // )).parse(input)
}

fn location_complement(input: &[u8]) -> IResult<&[u8], Location> {
    todo!("location_complement")
    // map(operator("complement", location), |l| {
    //     Location::Complement(Box::new(l))
    // }).parse(input)
}

fn location_gap(input: &[u8]) -> IResult<&[u8], Location> {
    todo!("location_gap")
    // let gap_length = alt((
    //     map(i64, GapLength::Known),
    //     value(GapLength::Unk100, tag("unk100")),
    //     value(GapLength::Unknown, tag("")),
    // ));
    // map(operator("gap", gap_length), Location::Gap).parse(input)
}

fn location_external(input: &[u8]) -> IResult<&[u8], Location> {
    let accession_parser = map_res(is_not(": \t\r\n"), str::from_utf8);
    let location = opt(preceded(
        tag(":"),
        alt((
            // everything except location_external
            location_span,
            location_compound,
            location_complement,
            location_between,
            location_single,
            location_gap,
        )),
    ));
    map(tuple((accession_parser, location)), |(accession, loc)| {
        Location::External(accession.to_owned(), loc.map(Box::new))
    }).parse(input)
}

pub fn location(input: &[u8]) -> IResult<&[u8], Location> {
    alt((
        location_span,
        location_compound,
        location_complement,
        location_between,
        location_single,
        location_gap,
        location_external,
    )).parse(input)
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_pos_single() {
        let s = "1\n";
        let p = location(s.as_bytes());
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
    fn test_complement_join() {
        let p = location(&b"join(complement(2..6),8)"[..]).unwrap();
        match p {
            (_, Location::Join(inner)) => match &*inner {
                [Location::Complement(inner), Location::Range((7, Before(false)), (8, After(false)))] => {
                    match inner.as_ref() {
                        Location::Range((1, Before(false)), (6, After(false))) => {}
                        x => {
                            panic!("{:?}", x);
                        }
                    }
                }
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
    fn test_range() {
        let s = "44..>579";
        assert_eq!(
            location(s.as_bytes()),
            Ok((
                (&b""[..]),
                Location::Range((43, Before(false)), (579, After(true)))
            ))
        )
    }

    #[test]
    fn test_pos_join() {
        let s = "join(1,3)";
        let p = location(s.as_bytes());
        let inner = vec![Location::simple_range(0, 1), Location::simple_range(2, 3)];
        assert_eq!(p, Ok(((&b""[..]), Location::Join(inner))));
    }
}
