use anyhow::{anyhow, Result};
use flate2::read::MultiGzDecoder;
use lazy_static::lazy_static;
use regex::bytes::Regex;
use std::str::FromStr;

use crate::helpers::{
    Address, Digram, Haplotype, LocalizedDigram, NodeId, Orientation, PathSegment, RawNodeId,
};
use std::{
    collections::{HashMap, HashSet},
    io::{self, BufRead, Read},
    path::PathBuf,
};

lazy_static! {
    static ref RE_WALK: Regex = Regex::new(r"([><])([!-;=?-~]+)").unwrap();
}

pub struct ByteLineReader<R: io::Read> {
    data: io::BufReader<R>,
}

impl<R: io::Read> ByteLineReader<R> {
    pub fn new(data: io::BufReader<R>) -> Self {
        Self { data }
    }
}

impl<R: io::Read> Iterator for ByteLineReader<R> {
    type Item = Vec<u8>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut buf = Vec::new();
        match self.data.read_until(b'\n', &mut buf) {
            Err(e) => {
                log::error!("Error while reading graph: {e:?} ");
                Some(buf)
            }
            Ok(1..) => Some(buf),
            Ok(0) => None,
        }
    }
}

fn bufreader_from_compressed(file: &PathBuf) -> Result<io::BufReader<Box<dyn Read>>> {
    let f = std::fs::File::open(file)?;
    if let Some(name) = file.file_name() {
        let reader: Box<dyn Read> = if name.to_str().unwrap().ends_with(".gz") {
            log::info!("assuming that {} is gzip compressed..", file.display());
            Box::new(MultiGzDecoder::new(f))
        } else {
            Box::new(f)
        };
        Ok(io::BufReader::new(reader))
    } else {
        Err(anyhow!("Filename {:?} is not a valid filename", file))
    }
}

pub fn parse_file_to_haplotypes(file: &PathBuf) -> Result<Vec<Haplotype>> {
    let node_ids_by_name = parse_node_ids(file, false)?;
    let data = bufreader_from_compressed(file)?;
    let line_reader = ByteLineReader::new(data);
    let haplotypes = line_reader
        .filter_map(|line| line_to_haplotype(line, &node_ids_by_name))
        .collect::<Vec<Haplotype>>();
    Ok(haplotypes)
}

fn line_to_haplotype(
    line: Vec<u8>,
    node_ids_by_name: &HashMap<Vec<u8>, RawNodeId>,
) -> Option<Haplotype> {
    if line.is_empty() {
        return None;
    }
    let first_char = line[0];
    match first_char {
        b'P' => path_to_haplotype(line, node_ids_by_name),
        b'W' => walk_to_haplotype(line, node_ids_by_name),
        _ => None,
    }
}

fn path_to_haplotype(
    line: Vec<u8>,
    node_ids_by_name: &HashMap<Vec<u8>, RawNodeId>,
) -> Option<Haplotype> {
    let (path_seg, buf_path_seg) = parse_path_identifier(&line);
    let haplotype = parse_path_seq(buf_path_seg, node_ids_by_name);
    Some(haplotype)
}

fn walk_to_haplotype(
    line: Vec<u8>,
    node_ids_by_name: &HashMap<Vec<u8>, RawNodeId>,
) -> Option<Haplotype> {
    let (path_seg, buf_path_seg) = parse_walk_identifier(&line);
    let haplotype = parse_walk_seq(buf_path_seg, node_ids_by_name);
    Some(haplotype)
}

fn parse_node_ids(gfa_file: &PathBuf, with_q: bool) -> Result<HashMap<Vec<u8>, RawNodeId>> {
    let mut node2id: HashMap<Vec<u8>, RawNodeId> = HashMap::default();

    log::info!("constructing indexes for node/edge IDs, node lengths, and P/W lines..");
    let mut node_id = 0; // important: id must be > 0, otherwise counting procedure will produce errors

    let mut buf = vec![];
    let mut data = bufreader_from_compressed(gfa_file)?;
    while data.read_until(b'\n', &mut buf).unwrap_or(0) > 0 {
        if buf[0] == b'S' || (with_q && buf[0] == b'Q') {
            let mut iter = buf[2..].iter();
            let offset = iter.position(|&x| x == b'\t').unwrap();
            if node2id
                .insert(buf[2..offset + 2].to_vec(), RawNodeId::new(node_id))
                .is_some()
            {
                panic!(
                    "Segment with ID {} occurs multiple times in GFA",
                    str::from_utf8(&buf[2..offset + 2]).unwrap()
                )
            }
            node_id += 1;
        }
        buf.clear();
    }

    log::info!("found: {} nodes", node2id.len());
    Ok(node2id)
}
fn parse_path_identifier(data: &[u8]) -> (PathSegment, &[u8]) {
    let mut iter = data.iter();

    let start = iter.position(|&x| x == b'\t').unwrap() + 1;
    let offset = iter.position(|&x| x == b'\t').unwrap();
    let path_name = str::from_utf8(&data[start..start + offset]).unwrap();
    (
        PathSegment::from_str(path_name),
        &data[start + offset + 1..],
    )
}

pub fn parse_walk_identifier(data: &[u8]) -> (PathSegment, &[u8]) {
    let mut six_col: Vec<&str> = Vec::with_capacity(6);

    let mut it = data.iter();
    let mut i = 0;
    for _ in 0..6 {
        let j = it.position(|x| x == &b'\t').unwrap();
        six_col.push(str::from_utf8(&data[i..i + j]).unwrap());
        i += j + 1;
    }

    let seq_start = match six_col[4] {
        "*" => None,
        a => Some(usize::from_str(a).unwrap()),
    };

    let seq_end = match six_col[5] {
        "*" => None,
        a => Some(usize::from_str(a).unwrap()),
    };

    let path_seg = PathSegment::new(
        six_col[1].to_string(),
        six_col[2].to_string(),
        six_col[3].to_string(),
        seq_start,
        seq_end,
    );

    (path_seg, &data[i..])
}

fn parse_path_seq(data: &[u8], node_ids_by_name: &HashMap<Vec<u8>, RawNodeId>) -> Haplotype {
    let mut haplotype = Vec::new();
    let mut prev_counter = 0;
    let mut nodes_visited_curr: HashSet<RawNodeId> = HashSet::new();
    let mut curr_counter = 0;
    let mut it = data.iter();
    let end = it
        .position(|x| x == &b'\t' || x == &b'\n' || x == &b'\r')
        .unwrap();

    let prev_node = data[..end]
        .split(|&x| x == b',')
        .take(1)
        .collect::<Vec<_>>()[0];
    let prev_node = prev_node.trim_ascii();
    let orientation = Orientation::from_path(prev_node[prev_node.len() - 1]);
    let prev_node = node_ids_by_name[&prev_node[..prev_node.len() - 1]];

    // Set the counter before setting the orientation to not take the orientation into account
    let mut prev_node = NodeId::new(prev_node, orientation);

    nodes_visited_curr.insert(prev_node.get_undirected());
    data[..end]
        .split(|&x| x == b',')
        .skip(1)
        .for_each(|current_node| {
            let current_node = current_node.trim_ascii();
            let orientation = Orientation::from_path(current_node[current_node.len() - 1]);
            let current_node = node_ids_by_name[&current_node[..current_node.len() - 1]];

            let current_node = NodeId::new(current_node, orientation);

            prev_counter = curr_counter;
            if nodes_visited_curr.contains(&current_node.get_undirected()) {
                curr_counter += 1;
                nodes_visited_curr.clear();
            }
            nodes_visited_curr.insert(current_node.get_undirected());

            let digram = Digram::new(prev_node, current_node);
            let address = Address::new(prev_counter, curr_counter);
            let local_digram = LocalizedDigram::new(digram, address);
            haplotype.push(local_digram);

            prev_node = current_node;
        });
    log::debug!("parsing path sequences of size {} bytes..", end);
    haplotype
}

fn parse_walk_seq(data: &[u8], node_ids_by_name: &HashMap<Vec<u8>, RawNodeId>) -> Haplotype {
    let mut haplotype = Vec::new();
    let mut prev_counter = 0;
    let mut nodes_visited_curr: HashSet<RawNodeId> = HashSet::new();
    let mut curr_counter = 0;
    let mut it = data.iter();
    let end = it
        .position(|x| x == &b'\t' || x == &b'\n' || x == &b'\r')
        .unwrap();

    let prev_node = RE_WALK
        .captures_iter(&data[..end])
        .take(1)
        .collect::<Vec<_>>();
    let orientation = Orientation::from_walk(prev_node[0][1][0]);
    let prev_node = node_ids_by_name[&prev_node[0][2]];

    let mut prev_node = NodeId::new(prev_node, orientation);

    nodes_visited_curr.insert(prev_node.get_undirected());

    RE_WALK.captures_iter(&data[..end]).skip(1).for_each(|m| {
        let orientation = Orientation::from_walk(m[1][0]);
        let current_node = node_ids_by_name[&m[2]];

        let current_node = NodeId::new(current_node, orientation);

        prev_counter = curr_counter;
        if nodes_visited_curr.contains(&current_node.get_undirected()) {
            curr_counter += 1;
            nodes_visited_curr.clear();
        }
        nodes_visited_curr.insert(current_node.get_undirected());

        let digram = Digram::new(prev_node, current_node);
        let address = Address::new(prev_counter, curr_counter);
        let local_digram = LocalizedDigram::new(digram, address);
        haplotype.push(local_digram);

        prev_node = current_node;
    });
    log::debug!("parsing walk sequences of size {} bytes..", end);
    haplotype
}
