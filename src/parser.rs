use core::str;
use flate2::read::MultiGzDecoder;
use indexmap::IndexMap;
use crate::path_segment::PathSegment;
use std::str::FromStr;
use crate::{NodeId, NeighborList, Digrams, ColorSet};

use lazy_static::lazy_static;
use priority_queue::PriorityQueue;
use regex::bytes::Regex;
use std::{
    collections::{BTreeSet, HashMap, HashSet},
    io::{BufRead, BufReader, Read},
};

lazy_static! {
    static ref RE_WALK: Regex = Regex::new(r"([><])([!-;=?-~]+)").unwrap();
}

pub fn bufreader_from_compressed_gfa(gfa_file: &str) -> BufReader<Box<dyn Read>> {
    log::info!("loading graph from {}", &gfa_file);
    let f = std::fs::File::open(gfa_file).expect("Error opening file");
    let reader: Box<dyn Read> = if gfa_file.ends_with(".gz") {
        log::info!("assuming that {} is gzip compressed..", &gfa_file);
        Box::new(MultiGzDecoder::new(f))
    } else {
        Box::new(f)
    };
    BufReader::new(reader)
}

pub fn canonize(x: NodeId, y: NodeId) -> (NodeId, NodeId) {
    if x > y || (x >> 1 == y >> 1 && x & 1 == 0) {
        (y ^ 1, x ^ 1)
    } else {
        (x, y)
    }
}

fn flip_digram(x: NodeId, y: NodeId) -> (NodeId, NodeId) {
    (y ^ 1, x ^ 1)
}

pub fn parse_gfa_paths_walks(
    gfa_file: &str,
    node_ids_by_name: &HashMap<Vec<u8>, NodeId>,
) -> (
    NeighborList,
    Digrams,
    HashMap<u64, PathSegment>,
) {
    log::info!("parsing path + walk sequences");
    let mut data = bufreader_from_compressed_gfa(gfa_file);

    let number_of_nodes = node_ids_by_name.len();
    let mut neighbors: Vec<BTreeSet<NodeId>> =
        vec![BTreeSet::new(); number_of_nodes * 2 + 1]; // Multiply by two to account for orientation
    let mut digrams: IndexMap<(NodeId, NodeId), ColorSet> = IndexMap::new();

    let mut path_id_to_path_segment: HashMap<u64, PathSegment> = HashMap::new();

    let mut buf = vec![];
    let mut path_id = 0;
    while data.read_until(b'\n', &mut buf).unwrap_or(0) > 0 {
        if buf[0] == b'P' || buf[0] == b'W' {
            let (path_seg, buf_path_seg) = match buf[0] {
                b'P' => parse_path_identifier(&buf),
                b'W' => parse_walk_identifier(&buf),
                _ => unreachable!(),
            };

            match buf[0] {
                b'P' => parse_path_seq(
                    buf_path_seg,
                    path_id,
                    node_ids_by_name,
                    &mut neighbors,
                    &mut digrams,
                ),
                b'W' => parse_walk_seq(
                    buf_path_seg,
                    path_id,
                    node_ids_by_name,
                    &mut neighbors,
                    &mut digrams,
                ),
                _ => unreachable!(),
            }
            path_id_to_path_segment.insert(path_id, path_seg);
            path_id += 1;
        }
        buf.clear();
    }
    let digrams: Vec<((NodeId, NodeId), ColorSet)> = digrams.into_iter().collect();
    let digrams: PriorityQueue<_, _> = PriorityQueue::from(digrams);
    (
        neighbors,
        digrams,
        path_id_to_path_segment,
    )
}

pub fn parse_path_seq(
    data: &[u8],
    path_id: u64,
    node_ids_by_name: &HashMap<Vec<u8>, NodeId>,
    neighbors: &mut [BTreeSet<NodeId>],
    digrams: &mut IndexMap<(NodeId, NodeId), ColorSet>,
) {
    log::debug!("Parsing path: {}", path_id);
    let mut nodes_visited: HashMap<NodeId, usize> = HashMap::new();
    let mut it = data.iter();
    let end = it
        .position(|x| x == &b'\t' || x == &b'\n' || x == &b'\r')
        .unwrap();

    let prev_node = data[..end]
        .split(|&x| x == b',')
        .take(1)
        .collect::<Vec<_>>()[0];
    let prev_node = prev_node.trim_ascii();
    let orientation = (prev_node[prev_node.len() - 1] == b'+') as NodeId;
    let prev_node = node_ids_by_name[&prev_node[..prev_node.len() - 1]];

    // Set the counter before setting the orientation to not take the orientation into account
    nodes_visited.insert(prev_node, 0);
    let mut prev_node = (prev_node << 1) ^ orientation;

    data[..end]
        .split(|&x| x == b',')
        .skip(1)
        .for_each(|current_node| {
            let current_node = current_node.trim_ascii();
            let orientation = (current_node[current_node.len() - 1] != b'+') as NodeId;
            let current_node = node_ids_by_name[&current_node[..current_node.len() - 1]];

            // Set the counter before setting the orientation to not take the orientation into account
            nodes_visited
                .entry(current_node)
                .and_modify(|counter| *counter += 1)
                .or_insert(0);
            let count = nodes_visited[&current_node];

            let current_node = (current_node << 1) ^ orientation;

            let (first_node, second_node) = canonize(prev_node, current_node);
            let (flipped_first_node, flipped_second_node) = flip_digram(first_node, second_node);

            neighbors[first_node as usize].insert(second_node);
            neighbors[flipped_first_node as usize].insert(flipped_second_node);

            let current_path_id = (count as u64) << 32 ^ path_id;

            digrams
                .entry((first_node, second_node))
                .and_modify(|c| {
                    c.0.insert(current_path_id);
                })
                .or_insert(ColorSet(HashSet::from([current_path_id])));

            prev_node = current_node;
        });

    log::debug!("parsing path sequences of size {} bytes..", end);
}

pub fn parse_walk_seq(
    data: &[u8],
    path_id: u64,
    node_ids_by_name: &HashMap<Vec<u8>, NodeId>,
    neighbors: &mut [BTreeSet<NodeId>],
    digrams: &mut IndexMap<(NodeId, NodeId), ColorSet>,
) {
    let mut nodes_visited: HashMap<NodeId, usize> = HashMap::new();
    let mut it = data.iter();
    let end = it
        .position(|x| x == &b'\t' || x == &b'\n' || x == &b'\r')
        .unwrap();

    let prev_node = RE_WALK
        .captures_iter(&data[..end])
        .take(1)
        .collect::<Vec<_>>();
    let orientation = (prev_node[0][1][0] == b'>') as NodeId;
    let prev_node = node_ids_by_name[&prev_node[0][2]];

    // Set the counter before setting the orientation to not take the orientation into account
    nodes_visited.insert(prev_node, 0);
    let mut prev_node = (prev_node << 1) | orientation;

    RE_WALK.captures_iter(&data[..end]).skip(1).for_each(|m| {
        let orientation = (m[1][0] != b'>') as NodeId;
        let current_node = node_ids_by_name[&m[2]];

        // Set the counter before setting the orientation to not take the orientation into account
        nodes_visited
            .entry(current_node)
            .and_modify(|counter| *counter += 1)
            .or_insert(0);
        let count = nodes_visited[&current_node];

        let current_node = (current_node << 1) ^ orientation;

        let (first_node, second_node) = canonize(prev_node, current_node);
        let (flipped_first_node, flipped_second_node) = flip_digram(first_node, second_node);

        neighbors[first_node as usize].insert(second_node);
        neighbors[flipped_first_node as usize].insert(flipped_second_node);

        let second_path_id = (count as u64) << 32 ^ path_id;

        digrams
            .entry((first_node, second_node))
            .and_modify(|c| {
                c.0.insert(second_path_id);
            })
            .or_insert(ColorSet(HashSet::from([second_path_id])));

        prev_node = current_node;
    });

    log::debug!("parsing path sequences of size {} bytes..", end);
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

pub fn parse_path_identifier(data: &[u8]) -> (PathSegment, &[u8]) {
    let mut iter = data.iter();

    let start = iter.position(|&x| x == b'\t').unwrap() + 1;
    let offset = iter.position(|&x| x == b'\t').unwrap();
    let path_name = str::from_utf8(&data[start..start + offset]).unwrap();
    (
        PathSegment::from_str(path_name),
        &data[start + offset + 1..],
    )
}

pub fn parse_node_ids(gfa_file: &str) -> HashMap<Vec<u8>, u64> {
    let mut node2id: HashMap<Vec<u8>, NodeId> = HashMap::default();

    log::info!("constructing indexes for node/edge IDs, node lengths, and P/W lines..");
    let mut node_id = 0; // important: id must be > 0, otherwise counting procedure will produce errors

    let mut buf = vec![];
    let mut data = bufreader_from_compressed_gfa(gfa_file);
    while data.read_until(b'\n', &mut buf).unwrap_or(0) > 0 {
        if buf[0] == b'S' {
            let mut iter = buf[2..].iter();
            let offset = iter.position(|&x| x == b'\t').unwrap();
            if node2id
                .insert(buf[2..offset + 2].to_vec(), node_id)
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
    node2id
}
