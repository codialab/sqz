use std::str::FromStr;
use std::{
    collections::{HashMap, HashSet},
    fmt,
};

use once_cell::sync::Lazy;
use regex::Regex;

static PATHID_PANSN: Lazy<Regex> =
    Lazy::new(|| Regex::new(r"^([^#]+)(#[^#]+)?(#[^#].*)?$").unwrap());
static PATHID_COORDS: Lazy<Regex> = Lazy::new(|| Regex::new(r"^(.+):([0-9]+)-([0-9]+)$").unwrap());

#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub struct AddressNumber(u32);

#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub struct Address(AddressNumber, AddressNumber);

impl Address {
    pub fn new(a: u32, b: u32) -> Self {
        Address(AddressNumber(a), AddressNumber(b))
    }

    pub fn flip(&self) -> Self {
        Address(self.1, self.0)
    }

    pub fn get_first(&self) -> AddressNumber {
        self.0
    }

    pub fn get_second(&self) -> AddressNumber {
        self.1
    }
}

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub enum Orientation {
    Forward,
    Backward,
}

impl Orientation {
    pub fn flip(&self) -> Self {
        match self {
            Orientation::Forward => Orientation::Backward,
            Orientation::Backward => Orientation::Forward,
        }
    }

    pub fn from_path(byte: u8) -> Self {
        match byte {
            b'+' => Orientation::Forward,
            b'-' => Orientation::Backward,
            _ => panic!(
                "Orientation {} ({}) is not valid for paths",
                byte, byte as char
            ),
        }
    }

    pub fn from_walk(byte: u8) -> Self {
        match byte {
            b'>' => Orientation::Forward,
            b'<' => Orientation::Backward,
            _ => panic!(
                "Orientation {} ({}) is not valid for walks",
                byte, byte as char
            ),
        }
    }
}

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub struct RawNodeId(u32);

impl RawNodeId {
    pub fn new(id: u32) -> Self {
        Self(id)
    }
}

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub struct NodeId(RawNodeId, Orientation);

impl NodeId {
    pub fn new(id_of_node: RawNodeId, orientation: Orientation) -> Self {
        Self(id_of_node, orientation)
    }

    pub fn flip(&self) -> Self {
        Self(self.0, self.1.flip())
    }

    pub fn same_node(&self, other: &Self) -> bool {
        self.0 == other.0
    }

    pub fn get_undirected(&self) -> RawNodeId {
        self.0
    }
}

#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub struct Digram(NodeId, NodeId);

impl Digram {
    pub fn new(u: NodeId, v: NodeId) -> Self {
        Self(u, v)
    }

    pub fn canonize(&self) -> Self {
        if self.0 > self.1 {
            Self(self.1.flip(), self.0.flip())
        } else {
            self.clone()
        }
    }

    pub fn is_canonical(&self) -> bool {
        if self.0 > self.1 {
            false
        } else {
            true
        }
    }
}

#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub struct CanonicalDigram(NodeId, NodeId);

impl CanonicalDigram {
    pub fn get_u(&self) -> NodeId {
        self.0
    }

    pub fn get_v(&self) -> NodeId {
        self.1
    }
}

impl From<Digram> for CanonicalDigram {
    fn from(item: Digram) -> Self {
        let canonical = item.canonize();
        Self(canonical.0, canonical.1)
    }
}

impl Into<Digram> for CanonicalDigram {
    fn into(self) -> Digram {
        Digram(self.0, self.1)
    }
}

#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct LocalizedDigram(Digram, Address);

impl LocalizedDigram {
    pub fn new(digram: Digram, address: Address) -> Self {
        Self(digram, address)
    }

    pub fn canonize(&self) -> Self {
        if self.0.is_canonical() {
            self.clone()
        } else {
            Self(self.0.canonize(), self.1.flip())
        }
    }

    pub fn is_canonical(&self) -> bool {
        if self.0.is_canonical() {
            true
        } else {
            false
        }
    }

    pub fn split_to_canonical(self) -> (CanonicalDigram, Address) {
        let canonical = self.canonize();
        (canonical.0.into(), canonical.1)
    }
}

pub type Haplotype = Vec<LocalizedDigram>;

pub struct D {
    inner: HashMap<CanonicalDigram, HashSet<(usize, Address)>>,
}

impl D {
    pub fn new() -> Self {
        Self {
            inner: HashMap::new(),
        }
    }

    pub fn insert(&mut self, index: &CanonicalDigram, value: (usize, Address)) {
        self.inner.entry(index.clone()).or_default().insert(value);
    }

    pub fn get_freq(&self) -> Freq {
        let mut freq = Freq::new();
        for (digram, occurrences) in &self.inner {
            freq.insert(digram.clone(), occurrences.len());
        }
        freq
    }
}

pub struct Freq {
    most_frequent_index: usize,
    inner: Vec<HashSet<CanonicalDigram>>,
}

impl Freq {
    pub fn new() -> Self {
        Self {
            most_frequent_index: 0,
            inner: Vec::new(),
        }
    }

    pub fn insert(&mut self, digram: CanonicalDigram, number_of_occurrences: usize) {
        if number_of_occurrences >= self.inner.len() {
            let new_len = number_of_occurrences + 1;
            self.inner.resize_with(new_len, Default::default);
            self.most_frequent_index = number_of_occurrences;
        }
        self.inner[number_of_occurrences].insert(digram);
    }

    pub fn get_most_freq(&mut self) -> CanonicalDigram {
        let elt = self.inner[self.most_frequent_index]
            .iter()
            .next()
            .cloned()
            .unwrap();
        let result = self.inner[self.most_frequent_index].take(&elt).unwrap();
        if self.inner[self.most_frequent_index].is_empty() {
            if self.most_frequent_index > 0 {
                self.most_frequent_index -= 1;
            }
        }
        result
    }
}

pub struct NeighborList {
    inner: HashMap<(NodeId, usize, AddressNumber), (NodeId, AddressNumber)>,
}

impl NeighborList {
    pub fn new() -> Self {
        Self {
            inner: HashMap::new(),
        }
    }

    pub fn insert(&mut self, key: (NodeId, usize, AddressNumber), value: (NodeId, AddressNumber)) {
        self.inner.insert(key, value);
    }
}

#[derive(Debug, Clone, PartialEq, PartialOrd, Hash, Eq, Ord)]
pub struct PathSegment {
    pub sample: String,
    pub haplotype: Option<String>,
    pub seqid: Option<String>,
    pub start: Option<usize>,
    pub end: Option<usize>,
}

impl PathSegment {
    pub fn new(
        sample: String,
        haplotype: String,
        seqid: String,
        start: Option<usize>,
        end: Option<usize>,
    ) -> Self {
        Self {
            sample,
            haplotype: Some(haplotype),
            seqid: Some(seqid),
            start,
            end,
        }
    }

    pub fn from_str(s: &str) -> Self {
        let mut res = PathSegment {
            sample: s.to_string(),
            haplotype: None,
            seqid: None,
            start: None,
            end: None,
        };

        if let Some(c) = PATHID_PANSN.captures(s) {
            let segments: Vec<&str> = c.iter().filter_map(|x| x.map(|y| y.as_str())).collect();
            // first capture group is the string itself
            match segments.len() {
                4 => {
                    res.sample = segments[1].to_string();
                    res.haplotype = Some(segments[2][1..].to_string());
                    match PATHID_COORDS.captures(&segments[3][1..]) {
                        None => {
                            res.seqid = Some(segments[3][1..].to_string());
                        }
                        Some(cc) => {
                            res.seqid = Some(cc.get(1).unwrap().as_str().to_string());
                            res.start = usize::from_str(cc.get(2).unwrap().as_str()).ok();
                            res.end = usize::from_str(cc.get(3).unwrap().as_str()).ok();
                            log::debug!("path has coordinates {} ", res);
                        }
                    }
                }
                3 => {
                    res.sample = segments[1].to_string();
                    match PATHID_COORDS.captures(&segments[2][1..]) {
                        None => {
                            res.haplotype = Some(segments[2][1..].to_string());
                        }
                        Some(cc) => {
                            res.haplotype = Some(cc.get(1).unwrap().as_str().to_string());
                            res.start = usize::from_str(cc.get(2).unwrap().as_str()).ok();
                            res.end = usize::from_str(cc.get(3).unwrap().as_str()).ok();
                            log::debug!("path has coordinates {} ", res);
                        }
                    }
                }
                2 => {
                    if let Some(cc) = PATHID_COORDS.captures(segments[1]) {
                        res.sample = cc.get(1).unwrap().as_str().to_string();
                        res.start = usize::from_str(cc.get(2).unwrap().as_str()).ok();
                        res.end = usize::from_str(cc.get(3).unwrap().as_str()).ok();
                        log::debug!("path has coordinates {}", res);
                    }
                }
                _ => (),
            }
        }
        res
    }

    #[allow(dead_code)]
    pub fn from_str_start_end(s: &str, start: usize, end: usize) -> Self {
        let mut segment = Self::from_str(s);
        segment.start = Some(start);
        segment.end = Some(end);
        segment
    }

    pub fn id(&self) -> String {
        if self.haplotype.is_some() {
            format!(
                "{}#{}{}",
                self.sample,
                self.haplotype.as_ref().unwrap(),
                if self.seqid.is_some() {
                    "#".to_owned() + self.seqid.as_ref().unwrap().as_str()
                } else {
                    "".to_string()
                }
            )
        } else if self.seqid.is_some() {
            format!(
                "{}#*#{}",
                self.sample,
                self.seqid.as_ref().unwrap().as_str()
            )
        } else {
            self.sample.clone()
        }
    }

    #[allow(dead_code)]
    pub fn clear_coords(&self) -> Self {
        Self {
            sample: self.sample.clone(),
            haplotype: self.haplotype.clone(),
            seqid: self.seqid.clone(),
            start: None,
            end: None,
        }
    }

    pub fn to_walk_string(&self) -> String {
        format!(
            "{}\t{}\t{}\t{}\t{}",
            self.sample,
            self.haplotype.as_ref().unwrap_or(&0.to_string()),
            self.seqid.as_ref().unwrap_or(&0.to_string()),
            self.start.map(|x| x.to_string()).unwrap_or("*".to_string()),
            self.end.map(|x| x.to_string()).unwrap_or("*".to_string())
        )
    }

    pub fn coords(&self) -> Option<(usize, usize)> {
        if self.start.is_some() && self.end.is_some() {
            Some((self.start.unwrap(), self.end.unwrap()))
        } else {
            None
        }
    }
}

impl fmt::Display for PathSegment {
    fn fmt(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
        if let Some((start, end)) = self.coords() {
            write!(formatter, "{}:{}-{}", self.id(), start, end)?;
        } else {
            write!(formatter, "{}", self.id())?;
        }
        Ok(())
    }
}
