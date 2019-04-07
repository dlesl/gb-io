use bio::alphabets::dna::revcomp;
use std::borrow::Cow;
use std::cmp;
use std::error::Error;
use std::fmt;
use std::io;
use std::io::Write;
use std::str;

use crate::errors::GbParserError;
use crate::reader::parse_location;
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
        if month >= 1 && month <= 12 && day >= 1 && day <= 31 {
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
    Unk100
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
/// nucleotide within the sequence. Ranges are inclusive. To specify a
/// range that wraps around, use join(x..last,1..y).
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(Debug, PartialEq, Clone)]
pub enum Location {
    /// Just a number
    Single(i64),
    /// n^n+1
    Between(i64, i64),
    /// [<]x..[>]y
    Span((i64, Before), (i64, After)),
    Complement(Box<Location>),
    Join(Vec<Location>),
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
    External(Location, Box<Error>),
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
    pub fn simple_span(a: i64, b: i64) -> Location {
        Location::Span((a, Before(false)), (b, After(false)))
    }

    /// Try to get the outermost bounds for a location. Returns the
    /// starting and finishing locations, as inclusive sequence
    /// coordinates.
    pub fn find_bounds(&self) -> Result<(i64, i64), LocationError> {
        match *self {
            Location::Span((a, _), (b, _)) => Ok((a, b)),
            Location::Complement(ref location) => location.find_bounds(),
            Location::Join(ref locations) => {
                let first = locations.first().ok_or(LocationError::Empty)?;
                let last = locations.last().unwrap();
                let (start, _) = first.find_bounds()?;
                let (_, end) = last.find_bounds()?;
                Ok((start, end))
            }
            Location::Single(p) => Ok((p, p)),
            Location::Between(a, b) => Ok((a, b)),
            Location::Order(ref locations)
            | Location::Bond(ref locations)
            | Location::OneOf(ref locations) => {
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
            Single(v) => Single(val(v)?),
            Between(a, b) => Between(val(a)?, val(b)?),
            Span((a, before), (b, after)) => Span((val(a)?, before), (val(b)?, after)),
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
            Single(a) => {
                if a >= start && a < end {
                    Some(self.clone())
                } else {
                    None
                }
            }
            Between(a, b) => {
                if (a >= start && a < end) || (b >= start && b < end) {
                    Some(self.clone())
                } else {
                    None
                }
            }
            Span((a, before), (b, after)) => {
                match (a >= start && a < end, b >= start && b < end) {
                    (true, true) => Some(simplify_shallow(Span((a, before), (b, after))).unwrap()),
                    (true, false) => {
                        Some(simplify_shallow(Span((a, before), (end - 1, After(false)))).unwrap())
                    } // TODO: this should be true?
                    (false, true) => {
                        Some(simplify_shallow(Span((start, Before(false)), (b, after))).unwrap())
                    } // TODO
                    (false, false) => {
                        // does the location span (start, end)?
                        if a <= start && b >= end - 1 {
                            Some(simplify_shallow(Location::simple_span(start, end - 1)).unwrap()) // TODO
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
            Location::Single(p) => format!("{}", p + 1),
            Location::Between(a, b) => format!("{}^{}", a + 1, b + 1),
            Location::Span((a, Before(before)), (b, After(after))) => {
                let before = if before { "<" } else { "" };
                let after = if after { ">" } else { "" };
                format!("{}{}..{}{}", before, a + 1, after, b + 1)
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
                Location::simple_span(start, end - 1)
            }
            Topology::Circular => {
                let (start, end) = self.unwrap_range(start, end);
                if end > self.len() {
                    assert!(end < self.len() * 2, "Range wraps around more than once!");
                    Location::Join(vec![
                        Location::simple_span(start, self.len() - 1),
                        Location::simple_span(0, end - self.len() - 1),
                    ])
                } else {
                    Location::simple_span(start, end - 1)
                }
            }
        };
        simplify(res).expect("invalid Location not possible here")
    }
    /// "Unwraps" a location on a circular sequence, so that coordinates that
    /// span the origin are replaced with locations beyond the origin.
    pub fn unwrap_location(&self, p: Location) -> Result<Location, LocationError> {
        let (first, last) = p.find_bounds()?;
        if first < 0 || last >= self.len() {
            return Err(LocationError::OutOfBounds(p));
        }
        if last < first && !self.is_circular() {
            return Err(LocationError::OutOfBounds(p));
        }
        let len = self.len();
        p.transform(&|p| Ok(p), &|v| {
            let res = if v < first { v + len } else { v };
            Ok(res)
        })
    }

    /// "Wraps" a location on a circular sequence, so that coordinates that
    /// extend beyond the end of the sequence are are wrapped to the origin.
    pub fn wrap_location(&self, p: Location) -> Result<Location, LocationError> {
        use Location::*;
        let res = p
            .transform(
                &|p| {
                    let res = match p {
                        Single(mut a) => {
                            while a >= self.len() {
                                a -= self.len();
                            }
                            Single(a)
                        }
                        Span((mut a, before), (mut b, after)) => {
                            while a >= self.len() {
                                a -= self.len();
                                b -= self.len();
                            }
                            if b < self.len() {
                                Span((a, before), (b, after))
                            } else {
                                Join(vec![
                                    Span((a, before), (self.len() - 1, After(false))),
                                    Span((0, Before(false)), (b - self.len(), after)),
                                ])
                            }
                        }
                        p => p,
                    };
                    Ok(res)
                },
                &|v| Ok(v),
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
            let moved = p.transform(&|p| Ok(p), &|v| Ok(v + shift)).unwrap(); // can't fail
            self.wrap_location(moved)
        } else {
            Ok(p.transform(&|p| Ok(p), &|v| Ok(v + shift)).unwrap())
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
                        Location::Span((a, Before(before)), (b, After(after))) => {
                            Location::Span((b, Before(after)), (a, After(before)))
                        }
                        Location::Between(a, b) => Location::Between(b, a),
                        p => p,
                    };
                    Ok(p)
                },
                &|v| Ok(self.len() - 1 - v),
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

    /// Extracts from `start` to `end`, keeping only features which
    /// fall entirely within this range. Note that `end` is not
    /// inclusive. Skips ambiguous features with a warning.
    pub fn extract_range_no_truncation(&self, start: i64, end: i64) -> Seq {
        debug!("Extract: {} to {}, len is {}", start, end, self.len());
        let (start, end) = self.unwrap_range(start, end);
        debug!("Normalised: {} to {}", start, end);
        let shift = -start;
        let mut features = Vec::new();
        for f in &self.features {
            if let Ok((x, y)) = f.location.find_bounds() {
                if (x < 0 || y < 0 || x > self.len() || y > self.len())
                    || (!self.is_circular() && y < x)
                {
                    warn!("Skipping feature with invalid location: {}", f.location);
                    continue;
                }
                let (mut x, mut y) = self.unwrap_range(x, y + 1); // to exclusive
                y -= 1; // back again
                if x < start {
                    x += self.len();
                    y += self.len();
                }
                if x >= start && y < end {
                    match self.relocate_feature(f.clone(), shift) {
                        Ok(f) => features.push(f),
                        Err(e) => warn!("Skipping feature with tricky location: {}", e),
                    }
                }
            }
        }
        Seq {
            features,
            seq: self.extract_range_seq(start, end).into(),
            ..Seq::empty()
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
        self.extract_location_with_fetcher(l, |_| Err(Box::new(NoFetcherError)))
    }

    /// Extract the sequence specified by `l`. If the location
    /// references an external sequence, `ext_fetcher` will be called
    /// with the name of this sequence to retrieve it.
    pub fn extract_location_with_fetcher<'a, F>(
        &self,
        l: &Location,
        mut ext_fetcher: F,
    ) -> Result<Vec<u8>, LocationError>
    where
        F: (FnMut(&str) -> Result<&'a Seq, Box<Error>>) + 'a,
    {
        self.extract_location_impl(l, &mut ext_fetcher)
    }

    fn extract_location_impl<'a, F>(
        &self,
        l: &Location,
        ext_fetcher: &mut F,
    ) -> Result<Vec<u8>, LocationError>
    where
        F: (FnMut(&str) -> Result<&'a Seq, Box<Error>>) + 'a,
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
                .get(usize_or(from)?..=usize_or(to)?)
                .ok_or_else(|| LocationError::OutOfBounds(l.clone()))?;
            Ok(s)
        };
        use Location::*;
        let res = match *l {
            Single(a) => get_range(a, a)?.into(),
            Span((a, _), (b, _)) => get_range(a, b)?.into(),
            Join(ref ls) => {
                let mut res = Vec::new();
                for l in ls {
                    res.extend_from_slice(&self.extract_location_impl(l, ext_fetcher)?);
                }
                res
            }
            Complement(ref l) => revcomp(self.extract_location_impl(l, ext_fetcher)?),
            External(ref name, ref ext_l) => {
                let ext_seq = ext_fetcher(name)
                    .map_err(|e| LocationError::External(l.clone(), e))?;
                if let Some(ext_l) = ext_l {
                    ext_seq.extract_location_impl(ext_l, ext_fetcher)?
                } else {
                    ext_seq.seq.clone()
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
                (Single(ref a), Single(b)) => {
                    if *a + 1 == b {
                        *last = Location::simple_span(*a, b);
                    } else if *a != b {
                        // ie. join(1,1) (can this happen?)
                        res.push(Single(b));
                    }
                }
                (Single(ref a), Span((c, Before(false)), d)) => {
                    if *a + 1 == c {
                        *last = Span((*a, Before(false)), d);
                    } else {
                        res.push(Span((c, Before(false)), d));
                    }
                }
                (Span(ref a, (ref b, After(false))), Single(d)) => {
                    if *b + 1 == d {
                        *last = Span(*a, (d, After(false)));
                    } else {
                        res.push(Single(d));
                    }
                }
                (Span(a, (ref b, After(false))), Span((c, Before(false)), d)) => {
                    if *b + 1 == c {
                        *last = Span(*a, d);
                    } else {
                        res.push(Span((c, Before(false)), d));
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
    p.transform(&simplify_shallow, &|v| Ok(v))
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
                    .all(|x| if let Complement(_) = x { true } else { false })
                {
                    xs = xs.into_iter().rev().map(Location::complement).collect();
                    Ok(Join(xs).complement())
                } else {
                    Ok(Join(xs))
                }
            }
        }
        Span((a, Before(false)), (b, After(false))) if a == b => Ok(Single(a)),
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
                Location::simple_span(0, 3),
                Location::simple_span(5, 7),
            ]),
            vec![Location::simple_span(0, 3), Location::simple_span(5, 7)]
        );
        assert_eq!(
            merge_adjacent(vec![
                Location::simple_span(0, 4),
                Location::simple_span(5, 7),
            ]),
            vec![Location::simple_span(0, 7)]
        );
        assert_eq!(
            merge_adjacent(vec![Location::simple_span(0, 4), Location::Single(5)]),
            vec![Location::simple_span(0, 5)]
        );
        assert_eq!(
            merge_adjacent(vec![Single(0), Location::simple_span(1, 5)]),
            vec![Location::simple_span(0, 5)]
        );
        assert_eq!(
            merge_adjacent(vec![Single(0), Single(1)]),
            vec![Location::simple_span(0, 1)]
        );
        assert_eq!(
            &Join(merge_adjacent(vec![
                Location::Span((1, Before(true)), (2, After(false))),
                Location::Span((3, Before(false)), (4, After(true)))
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
            s.relocate_location(Location::Single(5), -9).unwrap(),
            Location::Single(6)
        );
        let span1 = Location::simple_span(7, 9);
        let span2 = Location::simple_span(0, 3);
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
            &s.relocate_location(Location::simple_span(0, 9), 5)
                .unwrap()
                .to_gb_format(),
            "join(6..10,1..5)"
        );
        assert_eq!(
            &s.relocate_location(Location::simple_span(5, 7), 5)
                .unwrap()
                .to_gb_format(),
            "1..3"
        );
        assert_eq!(
            &s.relocate_location(
                Location::Join(vec![
                    Location::simple_span(7, 9),
                    Location::simple_span(0, 2)
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
                    Location::simple_span(8, 9),
                    Location::simple_span(0, 1)
                ]),
                2
            )
            .unwrap()
            .to_gb_format(),
            "1..4"
        );
        assert_eq!(
            &s.relocate_location(Location::simple_span(0, 2), 5)
                .unwrap()
                .to_gb_format(),
            "6..8"
        );
        assert_eq!(
            s.relocate_location(Location::Single(5), -9).unwrap(),
            Location::Single(6)
        );
        let span1 = Location::simple_span(7, 9);
        let span2 = Location::simple_span(0, 3);
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
    fn extract_seq_no_truncate() {
        use std::iter::repeat;
        let mut features = Vec::new();
        for i in 0..8 {
            let location = Location::simple_span(i, i + 2);
            features.push(Feature {
                // this _is_ inclusive
                location: location.clone(),
                kind: FeatureKind::from(format!("{}", location)),
                qualifiers: vec![],
            });
        }
        for i in 8..10 {
            features.push(Feature {
                location: Location::Join(vec![
                    Location::simple_span(i, 9),
                    Location::simple_span(0, i - 8),
                ]),
                kind: FeatureKind::from(format!("{}", i)),
                qualifiers: Vec::new(),
            });
        }
        let s = Seq {
            seq: repeat(b'A').take(10).collect(),
            topology: Topology::Linear,
            features,
            ..Seq::empty()
        };
        for i in 0..8 {
            // this _isn't_ inclusive
            let reg = s.extract_range_no_truncation(i, i + 3);
            assert_eq!(reg.features.len(), 1);
        }
        let s = Seq {
            topology: Topology::Circular,
            ..s
        };
        for i in -10..20 {
            let reg = s.extract_range_no_truncation(i, i + 3);
            assert_eq!(reg.features.len(), 1);
            println!(
                "{}, original location {}, current location {}",
                i, reg.features[0].kind, reg.features[0].location
            );
            assert_eq!(reg.features[0].location, Location::simple_span(0, 2));
        }
    }

    #[test]
    fn test_extract_linear() {
        let features = vec![Feature {
            location: Location::simple_span(0, 99),
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
                    simplify(Location::simple_span(0, j - 1)).unwrap()
                );
            }
        }
    }

    #[test]
    fn test_extract_circular() {
        let whole_seq = Feature {
            location: Location::simple_span(0, 9),
            kind: FeatureKind::from(""),
            qualifiers: Vec::new(),
        };
        let make_pos = |from: i64, to: i64| -> Location {
            if to >= 10 {
                Location::Join(vec![
                    simplify(Location::simple_span(from, 9)).unwrap(),
                    simplify(Location::simple_span(0, to - 10)).unwrap(),
                ])
            } else {
                simplify(Location::simple_span(from, to)).unwrap()
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
                            location: make_pos(i, i + j - 1),
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
                        simplify(Location::simple_span(0, j - 1)).unwrap()
                    );
                }
                println!("1 {:?}", res.features[0]);
                println!("{} {}", i, j);
                if i < 8 {
                    assert_eq!(
                        res.features[1].location,
                        simplify(Location::simple_span(0, j - 1)).unwrap()
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
                location: Location::simple_span(0, 3),
                kind: feature_kind!(""),
                qualifiers: vec![],
            }],
            ..Seq::empty()
        };
        let res = s.extract_range(4, 10);
        assert_eq!(res.features, vec![]);
        let res = s.extract_range(0, 1);
        assert_eq!(res.features.len(), 1);
        assert_eq!(&res.features[0].location, &Location::Single(0));
        let res = s.extract_range(3, 4);
        assert_eq!(res.features.len(), 1);
        assert_eq!(&res.features[0].location, &Location::Single(0));
        let res = s.extract_range(0, 10);
        assert_eq!(&res.features[0].location, &Location::simple_span(0, 3));
    }

    #[test]
    fn truncate() {
        assert_eq!(
            Location::Single(0).truncate(0, 1),
            Some(Location::Single(0))
        );
        assert_eq!(Location::Single(0).truncate(1, 2), None);
        assert_eq!(
            Location::simple_span(0, 2).truncate(1, 2),
            Some(Location::Single(1))
        );
        assert_eq!(
            Location::simple_span(0, 1).truncate(0, 2),
            Some(Location::simple_span(0, 1))
        );
        assert_eq!(
            Location::Complement(Box::new(Location::simple_span(0, 1))).truncate(0, 2),
            Some(Location::Complement(Box::new(Location::simple_span(0, 1))))
        );
        assert_eq!(
            Location::Complement(Box::new(Location::simple_span(0, 1))).truncate(10, 20),
            None
        );
        assert_eq!(Location::simple_span(0, 1).truncate(3, 4), None);
        let p = Location::Join(vec![
            Location::simple_span(0, 2),
            Location::simple_span(4, 6),
        ]);
        assert_eq!(p.truncate(0, 3), Some(Location::simple_span(0, 2)));
        assert_eq!(p.truncate(10, 30), None);
    }

    #[test]
    fn test_extract_circular_split() {
        let features = vec![Feature {
            location: Location::Join(vec![
                Location::simple_span(0, 2),
                Location::simple_span(4, 9),
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
                Location::simple_span(5, 9),
                Location::simple_span(0, 4),
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
                    location: Location::simple_span(2, 6),
                    kind: feature_kind!(""),
                    qualifiers: Vec::new(),
                },
                Feature {
                    location: Location::simple_span(0, 9),
                    kind: feature_kind!(""),
                    qualifiers: Vec::new(),
                },
                Feature {
                    location: Location::Join(vec![
                        Location::simple_span(7, 9),
                        Location::simple_span(0, 3),
                    ]),
                    kind: feature_kind!(""),
                    qualifiers: Vec::new(),
                },
                Feature {
                    location: Location::Single(0),
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
        assert_eq!(s.range_to_location(0, 10), Location::simple_span(0, 9));
        let s = Seq {
            topology: Topology::Circular,
            ..s
        };
        assert_eq!(s.range_to_location(5, 10), Location::simple_span(5, 9));
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
    fn unwrap_location() {
        let mut s = Seq {
            seq: "0123456789".into(),
            topology: Topology::Linear,
            features: Vec::new(),
            ..Seq::empty()
        };
        let location = Location::simple_span(2, 4);
        assert_eq!(
            s.unwrap_location(location).unwrap(),
            Location::simple_span(2, 4)
        );
        s.topology = Topology::Circular;
        let location = Location::simple_span(7, 3);
        assert_eq!(
            s.unwrap_location(location).unwrap(),
            Location::simple_span(7, 13)
        )
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
        assert_eq!(s.extract_location_with_fetcher(&l, |n| {
            assert_eq!(n, "TEST");
            Ok(&s_ext)
        }).unwrap(), vec![10, 11, 12, 13, 14]);
    }

    #[test]
    fn wrap_location() {
        let s = Seq {
            seq: (0..10).collect(),
            topology: Topology::Circular,
            ..Seq::empty()
        };
        assert_eq!(
            Location::Single(0),
            s.wrap_location(Location::Single(0)).unwrap()
        );
        assert_eq!(
            Location::Single(0),
            s.wrap_location(Location::Single(10)).unwrap()
        );
        assert_eq!(
            Location::Single(1),
            s.wrap_location(Location::Single(11)).unwrap()
        );
        assert_eq!(
            Location::simple_span(0, 1),
            s.wrap_location(Location::simple_span(10, 11)).unwrap()
        );
        assert_eq!(
            Location::Join(vec![Location::Single(9), Location::simple_span(0, 4)]),
            s.wrap_location(Location::simple_span(9, 14)).unwrap()
        );
        assert_eq!(
            &s.wrap_location(Location::Span((8, Before(false)), (10, After(true))))
                .unwrap()
                .to_gb_format(),
            "join(9..10,1..>1)"
        );
    }

    #[test]
    fn revcomp() {
        let make_seq = |locations: Vec<Location>| Seq {
            seq: "aaaaaaaaat".into(),
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
        assert_eq!(
            make_seq(vec![Location::simple_span(0, 1)])
                .revcomp()
                .features[0]
                .location,
            Location::Complement(Box::new(Location::simple_span(8, 9)))
        );
        assert_eq!(
            make_seq(vec![Location::Join(vec![
                Location::simple_span(0, 1),
                Location::simple_span(3, 4)
            ])])
            .revcomp()
            .features[0]
                .location,
            Location::Complement(Box::new(Location::Join(vec![
                Location::simple_span(5, 6),
                Location::simple_span(8, 9)
            ])))
        );
        assert_eq!(
            make_seq(vec![Location::Single(9)]).revcomp().features[0].location,
            Location::Complement(Box::new(Location::Single(0)))
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
