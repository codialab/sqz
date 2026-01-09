use anyhow::{anyhow, Result};
use std::hash::{BuildHasherDefault, DefaultHasher};
use std::mem;
use std::str::FromStr;
use std::{
    collections::{HashMap, HashSet},
    fmt,
};

use itertools::Itertools;
use once_cell::sync::Lazy;
use regex::Regex;

use crate::helpers::utils::{
    Address, AddressNumber, CanonicalDigram, LocalizedDigram, NodeId, Orientation, UndirectedNodeId,
};

pub mod digram_occurrences;
pub mod utils;

pub type DeterministicHashSet<T> = HashSet<T, BuildHasherDefault<DefaultHasher>>;
pub type DeterministicHashMap<T, U> = HashMap<T, U, BuildHasherDefault<DefaultHasher>>;

static PATHID_PANSN: Lazy<Regex> =
    Lazy::new(|| Regex::new(r"^([^#]+)(#[^#]+)?(#[^#].*)?$").unwrap());
static PATHID_COORDS: Lazy<Regex> = Lazy::new(|| Regex::new(r"^(.+):([0-9]+)-([0-9]+)$").unwrap());

#[derive(Debug)]
pub struct ReverseNodeRegistry {
    inner: DeterministicHashMap<UndirectedNodeId, Vec<u8>>,
}

impl ReverseNodeRegistry {
    pub fn get_directed_name(&self, node: NodeId) -> String {
        let o = match node.1 {
            Orientation::Forward => '>',
            Orientation::Backward => '<',
        };
        if let Some(node_name) = self.inner.get(&node.0) {
            return format!("{}{}", o, str::from_utf8(node_name).unwrap());
        } else {
            return format!("{}@{}", o, node.0 .0);
        }
    }

    pub fn get_name(&self, node: UndirectedNodeId) -> String {
        if let Some(node_name) = self.inner.get(&node) {
            return format!("{}", str::from_utf8(node_name).unwrap());
        } else {
            return format!("@{}", node.0);
        }
    }
}

#[derive(Debug)]
pub struct NodeRegistry {
    inner: DeterministicHashMap<Vec<u8>, UndirectedNodeId>,
    prefix: Vec<u8>,
    meta_node_number: u32,
}

impl NodeRegistry {
    pub fn new() -> Self {
        Self {
            inner: DeterministicHashMap::default(),
            prefix: vec![b'@'],
            meta_node_number: 1,
        }
    }

    #[allow(dead_code)]
    pub fn from(h: DeterministicHashMap<Vec<u8>, UndirectedNodeId>) -> Self {
        Self {
            inner: h,
            prefix: vec![b'@'],
            meta_node_number: 1,
        }
    }

    pub fn len(&self) -> usize {
        self.inner.len()
    }

    #[allow(dead_code)]
    pub fn with_prefix(prefix: Vec<u8>) -> Self {
        Self {
            inner: DeterministicHashMap::default(),
            prefix,
            meta_node_number: 1,
        }
    }

    pub fn insert(&mut self, name: Vec<u8>) -> Result<()> {
        self.get_inserted(name)
            .ok_or_else(|| anyhow!("Inserted node name already appeared"))?;
        Ok(())
    }

    fn get_inserted(&mut self, name: Vec<u8>) -> Option<UndirectedNodeId> {
        let new_node_id = self.inner.len();
        let new_node_id = UndirectedNodeId::new(new_node_id as u32);
        self.inner.insert(name, new_node_id);
        Some(new_node_id)
    }

    fn get_inserted_meta_node(&mut self, name: Vec<u8>) -> Option<UndirectedNodeId> {
        let new_node_id = self.inner.len();
        let new_node_id = UndirectedNodeId::new_meta_node(new_node_id as u32);
        self.inner.insert(name, new_node_id);
        Some(new_node_id)
    }

    #[allow(dead_code)]
    pub fn get_inserted_if_not_exists(&mut self, name: Vec<u8>) -> UndirectedNodeId {
        if self.inner.contains_key(&name) {
            self.get_id(&name)
        } else {
            self.get_inserted(name.clone())
                .expect(&format!("{:?}: {:?}", name, self.inner))
        }
    }

    pub fn get_id(&self, name: &[u8]) -> UndirectedNodeId {
        self.inner[name]
    }

    pub fn get_new_meta_node(&mut self) -> NodeId {
        let mut new_name = self.prefix.clone();
        new_name.extend(self.meta_node_number.to_string().into_bytes());
        let new_node_id = self
            .get_inserted_meta_node(new_name)
            .expect("Meta node name should be unique");
        self.meta_node_number += 1;
        let new_node_id = NodeId::new(new_node_id, Orientation::Forward);
        new_node_id
    }
}

impl Into<ReverseNodeRegistry> for NodeRegistry {
    fn into(self) -> ReverseNodeRegistry {
        ReverseNodeRegistry {
            inner: self.inner.into_iter().map(|(k, v)| (v, k)).collect(),
        }
    }
}

#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Occurrence(usize, Address);

impl fmt::Debug for Occurrence {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "[{}: {:?}]", self.0, self.1)
    }
}

impl Occurrence {
    pub fn new(index: usize, address: Address) -> Self {
        Self(index, address)
    }

    pub fn get_haplotype(&self) -> usize {
        self.0
    }

    pub fn get_address(&self) -> Address {
        self.1
    }

    #[allow(dead_code)]
    pub fn right_side_replacement(
        to_replace: HashSet<Occurrence, BuildHasherDefault<DefaultHasher>>,
        replacements: &HashMap<(usize, AddressNumber), AddressNumber>,
        flip: bool,
    ) -> HashSet<Occurrence, BuildHasherDefault<DefaultHasher>> {
        to_replace
            .into_iter()
            .map(|o| {
                if flip {
                    let replace = replacements[&(o.0, o.1 .1)];
                    Self(o.0, Address(o.1 .0, replace))
                } else {
                    let replace = replacements[&(o.0, o.1 .0)];
                    Self(o.0, Address(replace, o.1 .1))
                }
            })
            .collect()
    }

    #[allow(dead_code)]
    pub fn flip(&self) -> Self {
        Self(self.0, self.1.flip())
    }

    #[allow(dead_code)]
    pub fn flip_all(occurrences: &mut HashSet<Occurrence, BuildHasherDefault<DefaultHasher>>) {
        *occurrences = occurrences.iter().map(|o| o.flip()).collect();
    }

    pub fn split_self_loops(
        occurrences: HashSet<Occurrence, BuildHasherDefault<DefaultHasher>>,
    ) -> (Vec<Occurrence>, Vec<Occurrence>, Vec<Occurrence>) {
        let mut uv = Vec::new();
        let mut qq = Vec::new();
        let mut vq = Vec::new();
        let sections = Self::sectionize(occurrences);
        for mut section in sections {
            if section.len() % 2 == 0 {
                vq.push(
                    section
                        .pop()
                        .expect("Section needs to contain at least one element"),
                );
            }
            for mut chunk in &section.into_iter().chunks(2) {
                let first = chunk.next().expect("Chunk contains at least one element");
                uv.push(first.clone());
                if let Some(second) = chunk.next() {
                    println!("{:?} | {:?}", first, second);
                    let replace = Self(second.0, Address(first.1 .0, second.1 .1));
                    qq.push(replace);
                }
            }
        }
        (uv, qq, vq)
    }

    #[allow(dead_code)]
    pub fn split_pre_self_loop(
        mut vu: Vec<Occurrence>,
        uv: &HashSet<Occurrence, BuildHasherDefault<DefaultHasher>>,
    ) -> (
        Vec<Occurrence>,
        Vec<Occurrence>,
        Vec<Occurrence>,
        Vec<Occurrence>,
    ) {
        let mut qq = Vec::new();
        let mut vq = Vec::new();
        let mut qu = Vec::new();
        for occurrence_idx in 0..vu.len() {
            let mut matches_first_number = false;
            let mut matches_second_number = false;
            let occurrence = &vu[occurrence_idx];
            for uv_occurrence in uv {
                if occurrence.0 == uv_occurrence.0 {
                    if occurrence.1 .0 == uv_occurrence.1 .1 {
                        matches_first_number = true;
                    } else if occurrence.1 .1 == uv_occurrence.1 .0 {
                        matches_second_number = true;
                    }
                }
            }
            match (matches_first_number, matches_second_number) {
                (true, true) => qq.push(vu.swap_remove(occurrence_idx)),
                (false, true) => vq.push(vu.swap_remove(occurrence_idx)),
                (true, false) => qu.push(vu.swap_remove(occurrence_idx)),
                (false, false) => {}
            }
        }
        (qq, vq, qu, vu)
    }

    fn sectionize(
        occurrences: HashSet<Occurrence, BuildHasherDefault<DefaultHasher>>,
    ) -> Vec<Vec<Occurrence>> {
        let mut occurrences: Vec<Occurrence> = occurrences.into_iter().collect();
        occurrences.sort();
        let mut sections: Vec<Vec<Occurrence>> = Vec::new();
        let mut current_section = vec![occurrences[0].clone()];
        occurrences
            .iter()
            .tuple_windows()
            .for_each(|(first, second)| {
                if (first.0 == second.0) && (second.1 .0 == first.1 .1 || second.1 .1 == first.1 .0)
                {
                    current_section.push(second.clone());
                } else {
                    sections.push(mem::take(&mut current_section));
                    current_section = vec![second.clone()];
                }
            });
        sections.push(mem::take(&mut current_section));
        sections
    }
}

pub type Haplotype = (PathSegment, Vec<LocalizedDigram>);

#[derive(Clone, Debug)]
pub struct Freq {
    inner: Vec<DeterministicHashSet<CanonicalDigram>>,
}

impl Freq {
    pub fn new() -> Self {
        Self { inner: Vec::new() }
    }

    pub fn change_frequency(
        &mut self,
        digram: &CanonicalDigram,
        old_occurrence: usize,
        new_occurrence: usize,
    ) {
        self.remove(digram, old_occurrence);
        self.insert(digram.clone(), new_occurrence);
    }

    pub fn remove(&mut self, digram: &CanonicalDigram, number_of_occurrences: usize) {
        self.inner
            .get_mut(number_of_occurrences - 1)
            .expect("Freq needs to contain digram at occurrences to remove it")
            .remove(digram);
        // If the last one is empty, remove empty sets till a non-empty one is the last
        while self.inner.last().is_some_and(|set| set.is_empty()) {
            self.inner.pop();
        }
    }

    pub fn insert(&mut self, digram: CanonicalDigram, number_of_occurrences: usize) {
        if number_of_occurrences > self.inner.len() {
            self.inner
                .resize_with(number_of_occurrences, Default::default);
        }
        if number_of_occurrences == 0 {
            return;
        }
        self.inner[number_of_occurrences - 1].insert(digram);
    }

    pub fn get_most_freq(&self) -> Option<CanonicalDigram> {
        let elt = self.inner.last()?.iter().next().cloned();
        elt
    }
}

#[derive(Clone, Debug)]
pub struct NeighborList {
    inner: HashMap<(NodeId, usize, AddressNumber), (NodeId, AddressNumber)>,
}

impl NeighborList {
    pub fn new() -> Self {
        Self {
            inner: HashMap::new(),
        }
    }

    pub fn get_value(
        &self,
        node: NodeId,
        haplotype: usize,
        address_number: AddressNumber,
    ) -> Option<(NodeId, AddressNumber)> {
        self.inner.get(&(node, haplotype, address_number)).cloned()
    }

    pub fn insert(&mut self, key: (NodeId, usize, AddressNumber), value: (NodeId, AddressNumber)) {
        // Insert forward edge
        self.inner.insert(key.clone(), value.clone());

        // Insert reverse edge
        let rev_key = (value.0.flip(), key.1, value.1);
        let rev_value = (key.0.flip(), key.2);
        self.inner.insert(rev_key, rev_value);
    }

    pub fn remove(&mut self, key: (NodeId, usize, AddressNumber), value: (NodeId, AddressNumber)) {
        // Insert forward edge
        self.inner.remove(&key);

        // Insert reverse edge
        let rev_key = (value.0.flip(), key.1, value.1);
        self.inner.remove(&rev_key);
    }

    pub fn remove_occurrence(
        &mut self,
        digram: &CanonicalDigram,
        occurrence: &Occurrence,
        is_right_side: bool,
    ) {
        if is_right_side {
            let key = (digram.get_u(), occurrence.0, occurrence.1.get_first());
            let value = (digram.get_v(), occurrence.1.get_second());
            self.remove(key, value);
        } else {
            let key = (digram.get_v(), occurrence.0, occurrence.1.get_second());
            let value = (digram.get_u(), occurrence.1.get_first());
            self.remove(key, value);
        }
    }

    pub fn remove_occurrences(
        &mut self,
        digram: &CanonicalDigram,
        occurrences: &HashSet<Occurrence, BuildHasherDefault<DefaultHasher>>,
        is_right_side: bool,
    ) {
        for occurrence in occurrences {
            self.remove_occurrence(digram, occurrence, is_right_side);
        }
    }

    pub fn insert_occurrence(
        &mut self,
        digram: &CanonicalDigram,
        occurrence: Occurrence,
        is_right_side: bool,
    ) {
        if is_right_side {
            let key = (digram.get_u(), occurrence.0, occurrence.1.get_first());
            let value = (digram.get_v(), occurrence.1.get_second());
            self.insert(key, value);
        } else {
            let key = (digram.get_v(), occurrence.0, occurrence.1.get_second());
            let value = (digram.get_u(), occurrence.1.get_first());
            self.insert(key, value);
        }
    }

    pub fn insert_values(
        &mut self,
        digram: &CanonicalDigram,
        occurrences: Vec<Occurrence>,
        is_right_side: bool,
    ) {
        for occurrence in occurrences {
            self.insert_occurrence(digram, occurrence, is_right_side);
        }
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

    #[allow(dead_code)]
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

    #[allow(dead_code)]
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

    #[allow(dead_code)]
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
        write!(formatter, "{}", self.to_walk_string())?;
        Ok(())
    }
}
