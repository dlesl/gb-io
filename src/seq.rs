use std::borrow::{Borrow, Cow};
use std::cmp;
use std::error::Error;
use std::fmt;
use std::io;
use std::io::Write;
use std::str;

use crate::errors::GbParserError;
use crate::reader::parse_location;
use crate::dna::revcomp;
pub use crate::{FeatureKind, QualifierKey};

#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(Debug, PartialEq, Clone)]
/// A very simple Date struct so we don't have to rely on chrono
pub struct Date {
    year: i32,
    month: u32,
    day: u32,
}

#[derive(Debug, PartialEq, Clone)]
pub struct DateError;

impl Date {
    /// Construct from a calendar date, checks that the numbers look
    /// reasonable but nothing too exhaustive
    pub fn from_ymd(year: i32, month: u32, day: u32) -> Result<Date, DateError> {
        if (1..=12).contains(&month) && (1..=31).contains(&day) {
            Ok(Date { year, month, day })
        } else {
            Err(DateError)
        }
    }
    pub fn year(&self) -> i32 {
        self.year
    }
    pub fn month(&self) -> u32 {
        self.month
    }
    pub fn day(&self) -> u32 {
        self.day
    }
}

impl fmt::Display for Date {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let month = match self.month {
            1 => "JAN",
            2 => "FEB",
            3 => "MAR",
            4 => "APR",
            5 => "MAY",
            6 => "JUN",
            7 => "JUL",
            8 => "AUG",
            9 => "SEP",
            10 => "OCT",
            11 => "NOV",
            12 => "DEC",
            _ => unreachable!(),
        };
        write!(f, "{:02}-{}-{:04}", self.day, month, self.year)
    }
}

#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(Debug, PartialEq, Clone, Copy)]
pub enum GapLength {
    /// gap(n)
    Known(i64),
    /// gap()
    Unknown,
    /// gap(unk100)
    Unk100,
}

#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(Debug, PartialEq, Clone, Copy)]
pub struct Before(pub bool);
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(Debug, PartialEq, Clone, Copy)]
pub struct After(pub bool);

/// Represents a Genbank "location", used to specify the location of
/// features and in the CONTIG line. See
/// <http://www.insdc.org/files/feature_table.html> for a detailed
/// description of what they mean.
///
/// Note that locations specified here must always refer to a
/// nucleotide within the sequence. Ranges are inclusive in Genbank
/// format, but represented as exclusive ranges, using 0-based
/// indexing in this library. For example, `1..3` will be represented
/// as `Range((0, Before(false)), (4, After(false)))`.
///
/// To specify a range that wraps around on a circular sequence,
/// Genbank files use `join(x..last,1..y)`.
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(Debug, PartialEq, Clone)]
pub enum Location {
    /// Represents a range of positions, indicated with [<]x..[>]y in
    /// the Genbank file. If `<` or `>` are present, then `Before` or
    /// `After` will be true, respectively. This means that the
    /// feature starts before/extends beyond the end point. Genbank
    /// files represent locations consisting of a single position with
    /// a single number. In this library this is represented using
    /// `Range`, i.e. `1` becomes `Range((0, Before(false)), (1,
    /// After(false)))`.
    Range((i64, Before), (i64, After)),
    /// Represented as `n^n+1`: This means that the location is
    /// between the two _adjacent_ positions specified. On a circular
    /// sequence the last and first positions are also allowed.
    Between(i64, i64),
    /// INSDC: "Find the complement of the presented sequence in the
    /// span specified by "location" (i.e., read the complement of
    /// the presented strand in its 5'-to-3' direction)"
    Complement(Box<Location>),
    /// INSDC: "The indicated elements should be joined (placed
    /// end-to-end) to form one contiguous sequence"
    Join(Vec<Location>),
    /// INSDC: "The elements can be found in the specified order (5'
    /// to 3' direction), but nothing is implied about the
    /// reasonableness about joining them"
    Order(Vec<Location>),
    Bond(Vec<Location>),
    OneOf(Vec<Location>),
    External(String, Option<Box<Location>>),
    Gap(GapLength),
}

#[derive(Debug, Error)]
#[error(display = "Not configured to fetch external sequences")]
pub struct NoFetcherError;

#[derive(Debug, Error)]
pub enum LocationError {
    #[error(display = "Can't determine location due to ambiguity: {}", _0)]
    Ambiguous(Location),
    #[error(display = "Can't resolve external location `{}`: {}", _0, _1)]
    External(Location, Box<dyn Error>),
    #[error(display = "Recursion limit reached while processing: {}", _0)]
    Recursion(Location),
    // TODO: actually implement this
    #[error(display = "Empty location list encountered")]
    Empty,
    #[error(
        display = "Location refers to a location outside of the sequence: {}",
        _0
    )]
    OutOfBounds(Location),
}

impl Location {
    /// Convenience constructor for this commonly used variant
    pub fn simple_range(a: i64, b: i64) -> Location {
        Location::Range((a, Before(false)), (b, After(false)))
    }

    pub fn single(a: i64) -> Location {
        Location::simple_range(a, a + 1)
    }

    /// Try to get the start and end of a location. Returns the
    /// starting and finishing locations, as an exclusive range.
    pub fn find_bounds(&self) -> Result<(i64, i64), LocationError> {
        use Location::*;
        match *self {
            Range((a, _), (b, _)) => Ok((a, b)),
            Complement(ref location) => location.find_bounds(),
            Join(ref locations) | Order(ref locations) => {
                let first = locations.first().ok_or(LocationError::Empty)?;
                let last = locations.last().unwrap();
                let (start, _) = first.find_bounds()?;
                let (_, end) = last.find_bounds()?;
                Ok((start, end))
            }
            Between(a, b) => Ok((a, b + 1)),
            | Bond(ref locations)
            | OneOf(ref locations) => {
                let (left, right): (Vec<_>, Vec<_>) = locations
                    .iter()
                    .flat_map(Location::find_bounds) // This skips any Err values
                    .unzip();
                let min = left.into_iter().min().ok_or(LocationError::Empty)?;
                let max = right.into_iter().max().unwrap();
                Ok((min, max))
            }
            ref p => Err(LocationError::Ambiguous(p.clone())),
        }
    }

    // Only returns `Err` if one of the closures does
    fn transform<L, V>(self, loc: &L, val: &V) -> Result<Location, LocationError>
    where
        L: Fn(Location) -> Result<Location, LocationError>,
        V: Fn(i64) -> Result<i64, LocationError>,
    {
        loc(self)?.transform_impl(loc, val)
    }

    fn transform_impl<L, V>(self, loc: &L, val: &V) -> Result<Location, LocationError>
    where
        L: Fn(Location) -> Result<Location, LocationError>,
        V: Fn(i64) -> Result<i64, LocationError>,
    {
        use Location::*;
        let t_vec = |ps: Vec<Location>| {
            ps.into_iter()
                .map(|p| loc(p)?.transform_impl(loc, val))
                .collect::<Result<Vec<_>, _>>()
        };
        let res = match self {
            // Apply the location closure
            Complement(p) => Complement(Box::new(loc(*p)?.transform_impl(loc, val)?)),
            Order(ps) => Order(t_vec(ps)?),
            Bond(ps) => Bond(t_vec(ps)?),
            OneOf(ps) => OneOf(t_vec(ps)?),
            Join(ps) => Join(t_vec(ps)?),
            // Apply the value closure
            Between(a, b) => Between(val(a)?, val(b)?),
            Range((a, before), (b, after)) => Range((val(a)?, before), (val(b)?, after)),
            // We don't touch values here
            External(..) => self,
            Gap(..) => self,
        };
        Ok(res)
    }

    /// Truncates this location, limiting it to the given bounds.
    /// Note: `higher` is exclusive.
    /// `None` is returned if no part of the location lies within the bounds.
    pub fn truncate(&self, start: i64, end: i64) -> Option<Location> {
        use Location::*;
        let filter = |ps: &[Location]| {
            let res: Vec<_> = ps.iter().filter_map(|p| p.truncate(start, end)).collect();
            if res.is_empty() {
                None
            } else {
                Some(res)
            }
        };
        match *self {
            Between(a, b) => {
                if (a >= start && a < end) || (b >= start && b < end) {
                    Some(self.clone())
                } else {
                    None
                }
            }
            Range((a, before), (b, after)) => {
                match (a >= start && a < end, b > start && b <= end) {
                    (true, true) => Some(simplify_shallow(Range((a, before), (b, after))).unwrap()),
                    (true, false) => {
                        Some(simplify_shallow(Range((a, before), (end, After(false)))).unwrap())
                    } // TODO: this should be true?
                    (false, true) => {
                        Some(simplify_shallow(Range((start, Before(false)), (b, after))).unwrap())
                    } // TODO
                    (false, false) => {
                        // does the location span (start, end)?
                        if a <= start && b >= end {
                            Some(simplify_shallow(Location::simple_range(start, end)).unwrap()) // TODO
                        } else {
                            None
                        }
                    }
                }
            }
            Complement(ref a) => a.truncate(start, end).map(|a| Complement(Box::new(a))),
            Join(ref ps) => filter(ps).map(|v| simplify_shallow(Join(v)).unwrap()),
            OneOf(ref ps) => filter(ps).map(OneOf),
            Bond(ref ps) => filter(ps).map(Bond),
            Order(ref ps) => filter(ps).map(Order),
            External(_, _) | Gap(_) => Some(self.clone()),
        }
    }

    fn complement(self) -> Location {
        match self {
            Location::Complement(x) => *x,
            x => Location::Complement(Box::new(x)),
        }
    }

    pub fn to_gb_format(&self) -> String {
        fn location_list(locations: &[Location]) -> String {
            locations
                .iter()
                .map(Location::to_gb_format)
                .collect::<Vec<_>>()
                .join(",")
        }
        match *self {
            Location::Between(a, b) => format!("{}^{}", a + 1, b + 1),
            Location::Range((a, Before(before)), (b, After(after))) => {
                if !before && !after && b == a + 1 {
                    format!("{}", a + 1)
                } else {
                let before = if before { "<" } else { "" };
                let after = if after { ">" } else { "" };
                format!("{}{}..{}{}", before, a + 1, after, b)
                }
            }
            Location::Complement(ref p) => format!("complement({})", p.to_gb_format()),
            Location::Join(ref locations) => format!("join({})", location_list(locations)),
            Location::Order(ref locations) => format!("order({})", location_list(locations)),
            Location::Bond(ref locations) => format!("bond({})", location_list(locations)),
            Location::OneOf(ref locations) => format!("one-of({})", location_list(locations)),
            Location::External(ref name, Some(ref location)) => {
                format!("{}:{}", &name, location.to_gb_format())
            }
            Location::External(ref name, None) => format!("{}", name),
            Location::Gap(GapLength::Known(length)) => format!("gap({})", length),
            Location::Gap(GapLength::Unknown) => format!("gap()"),
            Location::Gap(GapLength::Unk100) => format!("gap(unk100)"),
        }
    }

    pub fn from_gb_format(s: &str) -> Result<Location, GbParserError> {
        parse_location(s.as_bytes())
    }
}

impl fmt::Display for Location {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.to_gb_format())
    }
}

/// An annotation, or "feature". Note that one qualifier key can occur multiple
/// times with distinct values. We store them in a `Vec` to preserve order. Some
/// qualifiers have no value.
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(Debug, PartialEq, Clone)]
pub struct Feature {
    pub kind: FeatureKind,
    pub location: Location,
    pub qualifiers: Vec<(QualifierKey, Option<String>)>,
}

impl Feature {
    /// Returns all the values for a given QualifierKey. Qualifiers with no
    /// value (ie. `/foo`) are ignored
    pub fn qualifier_values(&self, key: QualifierKey) -> impl Iterator<Item = &str> {
        self.qualifiers
            .iter()
            .filter(move |&(k, _)| k == &key)
            .filter_map(|&(_, ref v)| v.as_ref().map(String::as_str))
    }
}

#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(Debug, PartialEq, Clone)]
pub enum Topology {
    Linear,
    Circular,
}

impl fmt::Display for Topology {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let res = match *self {
            Topology::Linear => "linear",
            Topology::Circular => "circular",
        };
        write!(f, "{}", res)
    }
}

#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(Debug, PartialEq, Clone)]
pub struct Source {
    pub source: String,
    pub organism: Option<String>,
}

#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(Debug, PartialEq, Clone)]
pub struct Reference {
    pub description: String,
    pub authors: Option<String>,
    pub consortium: Option<String>,
    pub title: String,
    pub journal: Option<String>,
    pub pubmed: Option<String>,
    pub remark: Option<String>,
}

/// Maximum length for which the buffer holding the sequence will be
/// preallocated. This isn't a hard limit though... should it be?
#[doc(hidden)]
pub const REASONABLE_SEQ_LEN: usize = 500 * 1000 * 1000;

#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(Debug, PartialEq, Clone)]
pub struct Seq {
    /// Name as specified in the LOCUS line
    pub name: Option<String>,
    /// Whether this molecule is linear or circular
    pub topology: Topology,
    /// The date specified in the LOCUS line
    pub date: Option<Date>,
    /// Length as specified in the LOCUS line. Note that this may differ from
    /// `seq.len()`, especially if `contig` is `Some` and `seq` is empty. When
    /// parsing a file, if sequence data is provided we check that this value is
    /// equal to `seq.len()`
    pub len: Option<usize>,
    // TODO: should this be an option?
    pub molecule_type: Option<String>,
    pub division: String,
    pub definition: Option<String>,
    pub accession: Option<String>,
    pub version: Option<String>,
    pub source: Option<Source>,
    pub dblink: Option<String>,
    pub keywords: Option<String>,
    pub references: Vec<Reference>,
    pub comments: Vec<String>,
    #[cfg_attr(all(feature = "serde", feature = "serde_bytes"), serde(with = "serde_bytes"))]
    pub seq: Vec<u8>,
    pub contig: Option<Location>,
    pub features: Vec<Feature>,
}

impl Seq {
    /// Create a new, empty `Seq`
    pub fn empty() -> Seq {
        Seq {
            seq: vec![],
            contig: None,
            definition: None,
            accession: None,
            molecule_type: None,
            division: String::from("UNK"),
            dblink: None,
            keywords: None,
            source: None,
            version: None,
            comments: vec![],
            references: vec![],
            name: None,
            topology: Topology::Linear,
            date: None,
            len: None,
            features: vec![],
        }
    }

    pub fn is_circular(&self) -> bool {
        match self.topology {
            Topology::Circular => true,
            Topology::Linear => false,
        }
    }

    /// Returns the "actual" length of the sequence. Note that this may not be
    /// equal to self.seq.len(), in the following circumstances:
    /// - `self.seq` is empty and `self.contig` is not, and the corresponding
    /// file's LOCUS line specified a length. The returned value is i64 to
    /// simplifiy arithmetic with coords from `Location`
    pub fn len(&self) -> i64 {
        if let Some(len) = self.len {
            assert!(self.seq.is_empty() || len == self.seq.len());
            len as i64
        } else {
            self.seq.len() as i64
        }
    }

    /// "Normalises" an exclusive range on the chromosome, to make handling circular
    /// sequences simpler. For linear sequences, it doesn't do anything except
    /// panic unless `start` and `end` are both within the sequence, or if `end
    /// <= start`. For circular sequences, all values, including negative values,
    /// are allowed.
    ///
    /// If `end` <= `start`, the range is assumed to wrap around. The returned
    /// values will satisfy the following conditions:
    /// - `0 <= first < len`
    /// - `first < last`
    /// This means that in the case of a range which wraps around, `last` >= `len`.
    pub fn unwrap_range(&self, start: i64, end: i64) -> (i64, i64) {
        let len = self.len();
        match self.topology {
            Topology::Linear => {
                assert!(start <= end);
                assert!(start >= 0 && start < len && end >= 0 && end <= len);
                (start, end)
            }
            Topology::Circular => {
                let mut start = start;
                let mut end = end;

                if start >= end {
                    end += len;
                }

                while start >= len {
                    start -= len;
                    end -= len;
                }

                while start < 0 {
                    start += len;
                    end += len;
                }

                assert!(start >= 0 && start < len);
                assert!(end > start);

                (start, end)
            }
        }
    }
    /// Given a range on this sequence, returns a corresponding `Location`
    /// Note: this is a rust-style, exclusive range
    pub fn range_to_location(&self, start: i64, end: i64) -> Location {
        let res = match self.topology {
            Topology::Linear => {
                assert!(end > start);
                assert!(start < self.len() && end <= self.len());
                Location::simple_range(start, end)
            }
            Topology::Circular => {
                let (start, end) = self.unwrap_range(start, end);
                if end > self.len() {
                    assert!(end < self.len() * 2, "Range wraps around more than once!");
                    Location::Join(vec![
                        Location::simple_range(start, self.len()),
                        Location::simple_range(0, end - self.len()),
                    ])
                } else {
                    Location::simple_range(start, end)
                }
            }
        };
        simplify(res).expect("invalid Location not possible here")
    }

    /// "Wraps" a location on a circular sequence, so that coordinates that
    /// extend beyond the end of the sequence are are wrapped to the origin.
    pub fn wrap_location(&self, p: Location) -> Result<Location, LocationError> {
        use Location::*;
        let res = p
            .transform(
                &|p| {
                    let res = match p {
                        Range((mut a, before), (mut b, after)) => {
                            while a >= self.len() {
                                a -= self.len();
                                b -= self.len();
                            }
                            if b <= self.len() {
                                Range((a, before), (b, after))
                            } else {
                                Join(vec![
                                    Range((a, before), (self.len(), After(false))),
                                    Range((0, Before(false)), (b - self.len(), after)),
                                ])
                            }
                        }
                        p => p,
                    };
                    Ok(res)
                },
                &Ok,
            )
            .unwrap();
        simplify(res)
    }

    /// "Shifts" a feature forwards by `shift` NTs (can be negative)
    /// Note: If this fails you won't get the original `Feature`
    /// back. If this is important, you should clone first
    pub fn relocate_feature(&self, f: Feature, shift: i64) -> Result<Feature, LocationError> {
        let res = Feature {
            location: self.relocate_location(f.location, shift)?,
            ..f
        };
        Ok(res)
    }

    /// "Shifts" a location forwards by `shift` NTs (can be negative)
    /// Note: If this fails you won't get the original `Location`
    /// back. If this is important, you should clone first
    pub fn relocate_location(
        &self,
        p: Location,
        mut shift: i64,
    ) -> Result<Location, LocationError> {
        if self.is_circular() {
            while shift < 0 {
                shift += self.len();
            }
            while shift >= self.len() {
                shift -= self.len();
            }
            let moved = p.transform(&Ok, &|v| Ok(v + shift)).unwrap(); // can't fail
            self.wrap_location(moved)
        } else {
            Ok(p.transform(&Ok, &|v| Ok(v + shift)).unwrap())
        }
    }

    /// Used by `revcomp`
    fn revcomp_location(&self, p: Location) -> Result<Location, LocationError> {
        let p = p
            .transform(
                &|mut p| {
                    match p {
                        Location::Join(ref mut locations)
                        | Location::Order(ref mut locations)
                        | Location::Bond(ref mut locations)
                        | Location::OneOf(ref mut locations) => {
                            locations.reverse();
                        }
                        _ => (),
                    };
                    let p = match p {
                        Location::Range((a, Before(before)), (b, After(after))) => {
                            Location::Range((b, Before(after)), (a, After(before)))
                        }
                        Location::Between(a, b) => Location::Between(b + 1, a + 1),
                        p => p,
                    };
                    Ok(p)
                },
                &|v| Ok(self.len() - v),
            )
            .unwrap(); // can't fail
        simplify(p).map(Location::complement)
    }

    /// Note: If this fails you won't get the original `Feature`
    /// back. If this is important, you should clone first
    pub fn revcomp_feature(&self, f: Feature) -> Result<Feature, LocationError> {
        Ok(Feature {
            location: self.revcomp_location(f.location)?,
            ..f
        })
    }

    /// Returns the reverse complement of a `Seq`, skipping any features
    /// which can't be processed with a warning
    pub fn revcomp(&self) -> Seq {
        let mut features = Vec::with_capacity(self.features.len());
        for f in &self.features {
            match self.revcomp_feature(f.clone()) {
                Ok(f) => features.push(f),
                Err(e) => warn!("Encountered invalid feature location: {}", e),
            }
        }
        Seq {
            features,
            seq: revcomp(&self.seq),
            ..self.clone()
        }
    }

    /// Extracts just the sequence from `start` to `end`, taking into
    /// account circularity. Note that `end` is exclusive. Use
    /// this instead of `extract_range` if you don't need the
    /// features.
    pub fn extract_range_seq(&self, start: i64, end: i64) -> Cow<[u8]> {
        // Here we use usize everywhere for convenience, since we will never
        // have to deal with negative values
        let len = self.len() as usize;
        assert!(len == self.seq.len());
        let (start, end) = self.unwrap_range(start, end);
        let mut start = start as usize;
        let mut end = end as usize;
        if end <= len {
            Cow::from(&self.seq[start..end])
        } else {
            assert!(self.is_circular());
            let mut res = Vec::with_capacity(end - start);
            while end != 0 {
                let slice_end = cmp::min(end, len);
                end -= slice_end;
                res.extend(&self.seq[start..slice_end]);
                start = 0;
            }
            Cow::from(res)
        }
    }

    /// Extracts from `start` to `end`, truncating features that
    /// extend beyond this range.  Note that `end` is not
    /// inclusive. Skips ambiguous features with a warning.
    pub fn extract_range(&self, start: i64, end: i64) -> Seq {
        let (start, end) = self.unwrap_range(start, end);
        let mut shift = -start;
        if self.is_circular() {
            while shift < 0 {
                shift += self.len();
            }
            while shift > self.len() {
                shift -= self.len();
            }
        }
        let features = self
            .features
            .iter()
            .flat_map(
                |f| match self.relocate_location(f.location.clone(), shift) {
                    // let `truncate` filter locations outside the range
                    Ok(l) => l.truncate(0, end - start).map(|location| Feature {
                        location,
                        ..f.clone()
                    }),
                    Err(e) => {
                        warn!("Skipping feature, can't process invalid location: {}", e);
                        None
                    }
                },
            )
            .collect();
        Seq {
            features,
            seq: self.extract_range_seq(start, end).into(),
            ..Seq::empty()
        }
    }

    /// Extract the sequence specified by `l`. This version returns
    /// `Err(LocationError::External(_, NoFetcherError))` if it
    /// encounters a reference to an external sequence.
    ///
    /// See `extract_location_with_fetcher` for details
    pub fn extract_location(&self, l: &Location) -> Result<Vec<u8>, LocationError> {
        self.extract_location_with_fetcher(l, |_| Err::<Seq, _>(Box::new(NoFetcherError)))
    }

    /// Extract the sequence specified by `l`. If the location
    /// references an external sequence, `ext_fetcher` will be called
    /// with the name of this sequence to retrieve it. Since an external
    /// feature may be referenced multiple times, it might be best to
    /// return `&Seq` or similar from `ext_fetcher`.
    pub fn extract_location_with_fetcher<'a, F, S>(
        &self,
        l: &Location,
        mut ext_fetcher: F,
    ) -> Result<Vec<u8>, LocationError>
    where
        F: (FnMut(&str) -> Result<S, Box<dyn Error>>) + 'a,
        S: Borrow<Seq> + 'a,
    {
        self.extract_location_impl(l, &mut ext_fetcher)
    }

    fn extract_location_impl<'a, F, S>(
        &self,
        l: &Location,
        ext_fetcher: &mut F,
    ) -> Result<Vec<u8>, LocationError>
    where
        F: (FnMut(&str) -> Result<S, Box<dyn Error>>) + 'a,
        S: Borrow<Seq> + 'a,
    {
        let get_range = |from: i64, to: i64| -> Result<&[u8], LocationError> {
            let usize_or = |a: i64| -> Result<usize, LocationError> {
                if a < 0 {
                    Err(LocationError::OutOfBounds(l.clone()))
                } else {
                    Ok(a as usize)
                }
            };
            let s = self
                .seq
                .get(usize_or(from)?..usize_or(to)?)
                .ok_or_else(|| LocationError::OutOfBounds(l.clone()))?;
            Ok(s)
        };
        use Location::*;
        let res = match *l {
            Range((a, _), (b, _)) => get_range(a, b)?.into(),
            Join(ref ls) => {
                let mut res = Vec::new();
                for l in ls {
                    res.extend_from_slice(&self.extract_location_impl(l, ext_fetcher)?);
                }
                res
            }
            Complement(ref l) => revcomp(&self.extract_location_impl(l, ext_fetcher)?),
            External(ref name, ref ext_l) => {
                let ext_seq =
                    ext_fetcher(name).map_err(|e| LocationError::External(l.clone(), e))?;
                if let Some(ext_l) = ext_l {
                    ext_seq.borrow().extract_location_impl(ext_l, ext_fetcher)?
                } else {
                    ext_seq.borrow().seq.clone()
                }
            }
            _ => return Err(LocationError::Ambiguous(l.clone())),
        };
        Ok(res)
    }

    /// Returns a new `Seq`, rotated so that `origin` is at the start
    pub fn set_origin(&self, origin: i64) -> Seq {
        assert!(self.is_circular());
        assert!(origin < self.len());
        let rotated = self.extract_range_seq(origin, origin);
        Seq {
            seq: rotated.into(),
            features: self
                .features
                .iter()
                .cloned()
                .flat_map(|f| self.relocate_feature(f, -origin))
                .collect(),
            ..self.clone()
        }
    }

    pub fn write<T: Write>(&self, file: T) -> io::Result<()> {
        crate::writer::write(file, self)
    }
}

//TODO: should we merge adjacent locations when Before/After is set?
fn merge_adjacent(ps: Vec<Location>) -> Vec<Location> {
    use Location::*;
    let mut res: Vec<Location> = Vec::with_capacity(ps.len());
    for p in ps {
        if let Some(last) = res.last_mut() {
            match (&last, p) {
                (Range(a, (ref b, After(false))), Range((c, Before(false)), d)) => {
                    if *b == c {
                        *last = Range(*a, d);
                    } else {
                        res.push(Range((c, Before(false)), d));
                    }
                }
                (_, p) => res.push(p),
            }
        } else {
            res.push(p);
        }
    }
    res
}

fn flatten_join(v: Vec<Location>) -> Vec<Location> {
    let mut res = Vec::with_capacity(v.len());
    for p in v {
        match p {
            Location::Join(xs) => res.extend(flatten_join(xs)),
            p => res.push(p),
        }
    }
    res
}

/// This doesn't simplify everything yet...
/// TODO: return original Location somehow on failure
fn simplify(p: Location) -> Result<Location, LocationError> {
    p.transform(&simplify_shallow, &Ok)
}

fn simplify_shallow(p: Location) -> Result<Location, LocationError> {
    use Location::*;
    match p {
        Join(xs) => {
            if xs.is_empty() {
                return Err(LocationError::Empty);
            }
            let xs = flatten_join(xs);
            let mut xs = merge_adjacent(xs);
            assert!(!xs.is_empty());
            if xs.len() == 1 {
                // remove the join, we now have a new type of location
                // so we need to simplify again
                Ok(simplify_shallow(xs.pop().unwrap())?)
            } else {
                //if everything is 'complement', reverse the order and move it outside
                if xs
                    .iter()
                    .all(|x| matches!(x, Complement(_)))
                {
                    xs = xs.into_iter().rev().map(Location::complement).collect();
                    Ok(Join(xs).complement())
                } else {
                    Ok(Join(xs))
                }
            }
        }
        p => Ok(p),
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::tests::init;

    #[test]
    fn test_merge_adj() {
        use Location::*;
        assert_eq!(
            merge_adjacent(vec![
                Location::simple_range(0, 4),
                Location::simple_range(5, 8),
            ]),
            vec![Location::simple_range(0, 4), Location::simple_range(5, 8)]
        );
        assert_eq!(
            merge_adjacent(vec![
                Location::simple_range(0, 5),
                Location::simple_range(5, 8),
            ]),
            vec![Location::simple_range(0, 8)]
        );
        assert_eq!(
            merge_adjacent(vec![Location::simple_range(0, 5), Location::single(5)]),
            vec![Location::simple_range(0, 6)]
        );
        assert_eq!(
            merge_adjacent(vec![Location::single(0), Location::simple_range(1, 5)]),
            vec![Location::simple_range(0, 5)]
        );
        assert_eq!(
            merge_adjacent(vec![Location::single(0), Location::single(1)]),
            vec![Location::simple_range(0, 2)]
        );
        assert_eq!(
            &Join(merge_adjacent(vec![
                Location::Range((1, Before(true)), (3, After(false))),
                Location::Range((3, Before(false)), (5, After(true)))
            ]))
            .to_gb_format(),
            "join(<2..>5)"
        );
    }

    #[test]
    fn pos_relocate_circular() {
        let s = Seq {
            seq: b"0123456789".to_vec(),
            topology: Topology::Circular,
            ..Seq::empty()
        };
        assert_eq!(
            s.relocate_location(Location::single(5), -9).unwrap(),
            Location::single(6)
        );
        let span1 = Location::simple_range(7, 10);
        let span2 = Location::simple_range(0, 4);
        assert_eq!(span1.to_gb_format(), String::from("8..10"));
        let join = Location::Join(vec![span1, span2]);
        assert_eq!(
            s.relocate_location(join, -5).unwrap().to_gb_format(),
            String::from("3..9")
        );
    }

    #[test]
    fn relocate_location() {
        let s = Seq {
            seq: b"0123456789".to_vec(),
            topology: Topology::Circular,
            ..Seq::empty()
        };
        assert_eq!(
            &s.relocate_location(Location::simple_range(0, 10), 5)
                .unwrap()
                .to_gb_format(),
            "join(6..10,1..5)"
        );
        assert_eq!(
            &s.relocate_location(Location::simple_range(5, 8), 5)
                .unwrap()
                .to_gb_format(),
            "1..3"
        );
        assert_eq!(
            &s.relocate_location(
                Location::Join(vec![
                    Location::simple_range(7, 10),
                    Location::simple_range(0, 3)
                ]),
                2
            )
            .unwrap()
            .to_gb_format(),
            "join(10,1..5)"
        );
        assert_eq!(
            &s.relocate_location(
                Location::Join(vec![
                    Location::simple_range(8, 10),
                    Location::simple_range(0, 2)
                ]),
                2
            )
            .unwrap()
            .to_gb_format(),
            "1..4"
        );
        assert_eq!(
            &s.relocate_location(Location::simple_range(0, 3), 5)
                .unwrap()
                .to_gb_format(),
            "6..8"
        );
        assert_eq!(
            s.relocate_location(Location::single(5), -9).unwrap(),
            Location::single(6)
        );
        let span1 = Location::simple_range(7, 10);
        let span2 = Location::simple_range(0, 4);
        let join = Location::Join(vec![span1, span2]);
        assert_eq!(
            &s.relocate_location(join, -5).unwrap().to_gb_format(),
            "3..9"
        );
    }

    #[test]
    fn extract_range_seq() {
        let s = Seq {
            seq: b"0123456789".to_vec(),
            topology: Topology::Circular,
            ..Seq::empty()
        };
        assert_eq!(
            s.extract_range_seq(0, 10).into_owned(),
            b"0123456789".to_vec()
        );
        assert_eq!(
            s.extract_range_seq(0, 11).into_owned(),
            b"01234567890".to_vec()
        );
        assert_eq!(
            s.extract_range_seq(-1, 10).into_owned(),
            b"90123456789".to_vec()
        );
        assert_eq!(
            s.extract_range_seq(-1, 11).into_owned(),
            b"901234567890".to_vec()
        );
        assert_eq!(
            s.extract_range_seq(-1, 21).into_owned(),
            b"9012345678901234567890".to_vec()
        );
    }

    #[test]
    fn test_extract_linear() {
        let features = vec![Feature {
            location: Location::simple_range(0, 100),
            kind: FeatureKind::from(""),
            qualifiers: Vec::new(),
        }];
        let s = Seq {
            seq: ::std::iter::repeat(b'A').take(100).collect(),
            topology: Topology::Linear,
            features,
            ..Seq::empty()
        };
        for i in 0..91 {
            for j in 1..11 {
                let res = s.extract_range(i, i + j);
                assert_eq!(res.features.len(), 1);
                assert_eq!(
                    res.features[0].location,
                    simplify(Location::simple_range(0, j)).unwrap()
                );
            }
        }
    }

    #[test]
    fn test_extract_circular() {
        let whole_seq = Feature {
            location: Location::simple_range(0, 10),
            kind: FeatureKind::from(""),
            qualifiers: Vec::new(),
        };
        let make_pos = |from: i64, to: i64| -> Location {
            if to > 10 {
                Location::Join(vec![
                    simplify(Location::simple_range(from, 10)).unwrap(),
                    simplify(Location::simple_range(0, to - 10)).unwrap(),
                ])
            } else {
                simplify(Location::simple_range(from, to)).unwrap()
            }
        };

        for i in 0..10 {
            for j in 1..3 {
                let s = Seq {
                    seq: ::std::iter::repeat(b'A').take(10).collect(),
                    topology: Topology::Circular,
                    features: vec![
                        whole_seq.clone(),
                        Feature {
                            location: make_pos(i, i + j),
                            kind: FeatureKind::from(""),
                            qualifiers: Vec::new(),
                        },
                    ],
                    ..Seq::empty()
                };
                let res = s.extract_range(i, i + j);
                assert_eq!(res.features.len(), 2);
                if i < 8 {
                    assert_eq!(
                        res.features[0].location,
                        simplify(Location::simple_range(0, j)).unwrap()
                    );
                }
                println!("1 {:?}", res.features[0]);
                println!("{} {}", i, j);
                if i < 8 {
                    assert_eq!(
                        res.features[1].location,
                        simplify(Location::simple_range(0, j)).unwrap()
                    );
                }
                println!("2 {:?}", res.features[1]);
            }
        }
    }

    #[test]
    fn extract_exclude_features() {
        let s = Seq {
            topology: Topology::Circular,
            seq: (0..10).collect(),
            features: vec![Feature {
                location: Location::simple_range(0, 4),
                kind: feature_kind!(""),
                qualifiers: vec![],
            }],
            ..Seq::empty()
        };
        let res = s.extract_range(4, 10);
        assert_eq!(res.features, vec![]);
        let res = s.extract_range(0, 1);
        assert_eq!(res.features.len(), 1);
        assert_eq!(&res.features[0].location, &Location::single(0));
        let res = s.extract_range(3, 4);
        assert_eq!(res.features.len(), 1);
        assert_eq!(&res.features[0].location, &Location::single(0));
        let res = s.extract_range(0, 10);
        assert_eq!(&res.features[0].location, &Location::simple_range(0, 4));
    }

    #[test]
    fn truncate() {
        assert_eq!(
            Location::single(0).truncate(0, 1),
            Some(Location::single(0))
        );
        assert_eq!(Location::single(0).truncate(1, 2), None);
        assert_eq!(
            Location::simple_range(0, 3).truncate(1, 2),
            Some(Location::single(1))
        );
        assert_eq!(
            Location::simple_range(0, 2).truncate(0, 2),
            Some(Location::simple_range(0, 2))
        );
        assert_eq!(
            Location::Complement(Box::new(Location::simple_range(0, 2))).truncate(0, 2),
            Some(Location::Complement(Box::new(Location::simple_range(0, 2))))
        );
        assert_eq!(
            Location::Complement(Box::new(Location::simple_range(0, 2))).truncate(10, 20),
            None
        );
        assert_eq!(Location::simple_range(0, 2).truncate(3, 4), None);
        let p = Location::Join(vec![
            Location::simple_range(0, 3),
            Location::simple_range(4, 7),
        ]);
        assert_eq!(p.truncate(0, 3), Some(Location::simple_range(0, 3)));
        assert_eq!(p.truncate(10, 30), None);
    }

    #[test]
    fn test_extract_circular_split() {
        let features = vec![Feature {
            location: Location::Join(vec![
                Location::simple_range(0, 2),
                Location::simple_range(4, 9),
            ]),
            kind: FeatureKind::from(""),
            qualifiers: Vec::new(),
        }];
        let s = Seq {
            seq: (0..10).collect(),
            topology: Topology::Circular,
            features,
            ..Seq::empty()
        };

        for i in 0..10 {
            let res = s.extract_range(i, i + 10);
            println!("{:?}", res.features);
        }
    }

    #[test]
    fn extract_range_wrapped_linear() {
        let features = vec![Feature {
            location: Location::Join(vec![
                Location::simple_range(5, 9),
                Location::simple_range(0, 4),
            ]),
            kind: FeatureKind::from(""),
            qualifiers: Vec::new(),
        }];
        let s = Seq {
            seq: (0..10).collect(),
            features,
            ..Seq::empty()
        };
        let res = s.extract_range(2, 8);
        assert_eq!(res.features.len(), 1);
        println!("{:?}", res.features[0].location);
    }
    #[test]
    fn set_origin() {
        init();
        let seq = Seq {
            seq: "0123456789".into(),
            features: vec![
                Feature {
                    location: Location::simple_range(2, 7),
                    kind: feature_kind!(""),
                    qualifiers: Vec::new(),
                },
                Feature {
                    location: Location::simple_range(0, 10),
                    kind: feature_kind!(""),
                    qualifiers: Vec::new(),
                },
                Feature {
                    location: Location::Join(vec![
                        Location::simple_range(7, 10),
                        Location::simple_range(0, 4),
                    ]),
                    kind: feature_kind!(""),
                    qualifiers: Vec::new(),
                },
                Feature {
                    location: Location::single(0),
                    kind: feature_kind!(""),
                    qualifiers: Vec::new(),
                },
            ],
            topology: Topology::Circular,
            ..Seq::empty()
        };
        for i in 1..9 {
            println!("******** {}", i);
            let rotated = seq.set_origin(i);
            println!("rotated {:?}", rotated.features);
            let rotated2 = rotated.set_origin(10 - i);
            println!("rotated2: {:?}", rotated2.features);
            assert_eq!(rotated2.features, seq.features);
        }
    }
    #[test]
    fn range_to_location() {
        let s = Seq {
            seq: "0123456789".into(),
            topology: Topology::Linear,
            ..Seq::empty()
        };
        assert_eq!(s.range_to_location(0, 10), Location::simple_range(0, 10));
        let s = Seq {
            topology: Topology::Circular,
            ..s
        };
        assert_eq!(s.range_to_location(5, 10), Location::simple_range(5, 10));
        assert_eq!(s.range_to_location(5, 11).to_gb_format(), "join(6..10,1)");
        assert_eq!(
            s.range_to_location(5, 15).to_gb_format(),
            "join(6..10,1..5)"
        );
    }
    #[test]
    fn unwrap_range_linear() {
        let s = Seq {
            seq: "01".into(),
            ..Seq::empty()
        };
        assert_eq!(s.unwrap_range(0, 1), (0, 1));
        assert_eq!(s.unwrap_range(0, 2), (0, 2));
    }

    #[test]
    #[should_panic]
    fn unwrap_range_linear_bad() {
        let s = Seq {
            seq: "01".into(),
            ..Seq::empty()
        };
        let _ = s.unwrap_range(0, 3);
    }

    #[test]
    fn unwrap_range_circular() {
        let s = Seq {
            seq: "01".into(),
            topology: Topology::Circular,
            ..Seq::empty()
        };
        assert_eq!(s.unwrap_range(0, 1), (0, 1));
        assert_eq!(s.unwrap_range(0, 2), (0, 2));
        assert_eq!(s.unwrap_range(0, 3), (0, 3));
        assert_eq!(s.unwrap_range(-1, 0), (1, 2));
        assert_eq!(s.unwrap_range(1, 2), (1, 2));
        assert_eq!(s.unwrap_range(2, 3), (0, 1));
        assert_eq!(s.unwrap_range(-2, -1), (0, 1));
    }

    #[test]
    fn extract_location() {
        let s = Seq {
            seq: (0..10).collect(),
            topology: Topology::Linear,
            ..Seq::empty()
        };
        let p = |l| Location::from_gb_format(l).unwrap();
        let e = |l| s.extract_location(&p(l)).unwrap();
        assert_eq!(e("1"), vec![0]);
        assert_eq!(e("2..6"), vec![1, 2, 3, 4, 5]);
        assert_eq!(e("join(complement(2..6),8)"), vec![5, 4, 3, 2, 1, 7]);
        assert_eq!(e("1..2"), vec![0, 1]);
        assert_eq!(e("complement(1..2)"), vec![1, 0]);
        assert_eq!(e("complement(join(8..10,1..2))"), vec![1, 0, 9, 8, 7]);
        assert_eq!(
            e("join(complement(1..2),complement(8..10))"),
            vec![1, 0, 9, 8, 7]
        );
    }

    #[test]
    fn extract_location_external() {
        let s = Seq {
            seq: (0..10).collect(),
            ..Seq::empty()
        };
        let s_ext = Seq {
            seq: (10..20).collect(),
            ..Seq::empty()
        };
        let l = Location::from_gb_format("TEST:1..5").unwrap();
        assert_eq!(
            s.extract_location_with_fetcher(&l, |n| {
                assert_eq!(n, "TEST");
                Ok(&s_ext)
            })
            .unwrap(),
            vec![10, 11, 12, 13, 14]
        );
    }

    #[test]
    fn wrap_location() {
        let s = Seq {
            seq: (0..10).collect(),
            topology: Topology::Circular,
            ..Seq::empty()
        };
        assert_eq!(
            Location::single(0),
            s.wrap_location(Location::single(0)).unwrap()
        );
        assert_eq!(
            Location::single(0),
            s.wrap_location(Location::single(10)).unwrap()
        );
        assert_eq!(
            Location::single(1),
            s.wrap_location(Location::single(11)).unwrap()
        );
        assert_eq!(
            Location::simple_range(0, 2),
            s.wrap_location(Location::simple_range(10, 12)).unwrap()
        );
        assert_eq!(
            Location::Join(vec![Location::single(9), Location::simple_range(0, 5)]),
            s.wrap_location(Location::simple_range(9, 15)).unwrap()
        );
        assert_eq!(
            &s.wrap_location(Location::Range((8, Before(false)), (11, After(true))))
                .unwrap()
                .to_gb_format(),
            "join(9..10,1..>1)"
        );
    }

    #[test]
    fn revcomp() {
        let make_seq = |locations: Vec<Location>| Seq {
            seq: (0..10).collect(),
            topology: Topology::Linear,
            features: locations
                .into_iter()
                .map(|p| Feature {
                    location: p,
                    kind: feature_kind!(""),
                    qualifiers: Vec::new(),
                })
                .collect(),
            ..Seq::empty()
        };
        let s = make_seq(vec![]);
        assert_eq!(
            s.revcomp_location(Location::single(0)).unwrap(),
            Location::Complement(Box::new(Location::single(9)))
        );
        assert_eq!(
            s.revcomp_location(Location::Between(0, 1)).unwrap(),
            Location::Complement(Box::new(Location::Between(8,9)))
        );
        assert_eq!(
            make_seq(vec![Location::simple_range(0, 2)])
                .revcomp()
                .features[0]
                .location,
            Location::Complement(Box::new(Location::simple_range(8, 10)))
        );
        assert_eq!(
            make_seq(vec![Location::Join(vec![
                Location::simple_range(0, 2),
                Location::simple_range(3, 5)
            ])])
            .revcomp()
            .features[0]
                .location,
            Location::Complement(Box::new(Location::Join(vec![
                Location::simple_range(5, 7),
                Location::simple_range(8, 10)
            ])))
        );
        assert_eq!(
            make_seq(vec![Location::single(9)]).revcomp().features[0].location,
            Location::Complement(Box::new(Location::single(0)))
        );
    }
    #[test]
    fn test_simplify() {
        let s = |a, b| {
            assert_eq!(
                simplify(parse_location(a).unwrap()).unwrap().to_gb_format(),
                b
            );
        };
        s(b"1..5", "1..5");
        s(b"join(1..2,2..5)", "join(1..2,2..5)");
        s(b"join(1..2,3..5)", "1..5");
        s(b"join(join(1..2,3..4),5)", "1..5");
        s(
            b"join(complement(1..2),complement(4..5))",
            "complement(join(4..5,1..2))",
        );
    }
}
