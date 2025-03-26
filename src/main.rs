mod path_segment;

use clap::Parser;
use core::str;
use flate2::read::MultiGzDecoder;
use indexmap::IndexMap;
use lazy_static::lazy_static;
use path_segment::PathSegment;
use priority_queue::PriorityQueue;
use regex::bytes::Regex;
use std::str::FromStr;
use std::{
    collections::{BTreeSet, HashMap, HashSet},
    io::{BufRead, BufReader, Read},
};

const MAX_OCCURENCES: usize = 2;

type NodeId = u64;
type NeighborList = Vec<BTreeSet<NodeId>>;
type Digrams = PriorityQueue<(NodeId, NodeId), ColorSet>;
type Rules = HashMap<NodeId, (NodeId, NodeId)>;

lazy_static! {
    static ref RE_WALK: Regex = Regex::new(r"([><])([!-;=?-~]+)").unwrap();
}

#[derive(PartialEq, Eq, Debug, Clone)]
pub struct ColorSet(HashSet<u64>);

impl ColorSet {
    fn intersection(&self, other: &ColorSet) -> ColorSet {
        ColorSet(
            self.0
                .intersection(&other.0)
                .copied()
                .collect::<HashSet<_>>(),
        )
    }

    fn difference(&self, other: &ColorSet) -> ColorSet {
        ColorSet(self.0.difference(&other.0).copied().collect::<HashSet<_>>())
    }
}

impl PartialOrd for ColorSet {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for ColorSet {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.0.len().cmp(&other.0.len())
    }
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

fn canonize_digram(x: NodeId, y: NodeId) -> (NodeId, NodeId) {
    if x > y || (x >> 1 == y >> 1 && x & 1 == 0) {
        (y ^ 1, x ^ 1)
    } else {
        (x, y)
    }
}

pub fn parse_gfa_paths_walks(
    gfa_file: &str,
    node_ids_by_name: &HashMap<Vec<u8>, NodeId>,
) -> (
    NeighborList,
    NeighborList,
    Digrams,
    HashMap<u64, PathSegment>,
) {
    log::info!("parsing path + walk sequences");
    let mut data = bufreader_from_compressed_gfa(gfa_file);

    let number_of_nodes = node_ids_by_name.len();
    let mut forward_neighbors: Vec<BTreeSet<NodeId>> =
        vec![BTreeSet::new(); number_of_nodes * 2 + 1]; // Multiply by two to account for orientation
    let mut backward_neighbors: Vec<BTreeSet<NodeId>> =
        vec![BTreeSet::new(); number_of_nodes * 2 + 1];
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
                    &mut forward_neighbors,
                    &mut backward_neighbors,
                    &mut digrams,
                ),
                b'W' => parse_walk_seq(
                    buf_path_seg,
                    path_id,
                    node_ids_by_name,
                    &mut forward_neighbors,
                    &mut backward_neighbors,
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
        forward_neighbors,
        backward_neighbors,
        digrams,
        path_id_to_path_segment,
    )
}

pub fn parse_path_seq(
    data: &[u8],
    path_id: u64,
    node_ids_by_name: &HashMap<Vec<u8>, NodeId>,
    forward_neighbors: &mut [BTreeSet<NodeId>],
    backward_neighbors: &mut [BTreeSet<NodeId>],
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
            let orientation = (current_node[current_node.len() - 1] == b'+') as NodeId;
            let current_node = node_ids_by_name[&current_node[..current_node.len() - 1]];

            // Set the counter before setting the orientation to not take the orientation into account
            nodes_visited
                .entry(current_node)
                .and_modify(|counter| *counter += 1)
                .or_insert(0);
            let count = nodes_visited[&current_node];

            let current_node = (current_node << 1) ^ orientation;

            let (first_node, second_node) = canonize_digram(prev_node, current_node);

            forward_neighbors[first_node as usize].insert(second_node);
            backward_neighbors[second_node as usize].insert(first_node);

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
    forward_neighbors: &mut [BTreeSet<NodeId>],
    backward_neighbors: &mut [BTreeSet<NodeId>],
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
        let orientation = (m[1][0] == b'>') as NodeId;
        let current_node = node_ids_by_name[&m[2]];

        // Set the counter before setting the orientation to not take the orientation into account
        nodes_visited
            .entry(current_node)
            .and_modify(|counter| *counter += 1)
            .or_insert(0);
        let count = nodes_visited[&current_node];

        let current_node = (current_node << 1) ^ orientation;

        let (first_node, second_node) = canonize_digram(prev_node, current_node);

        forward_neighbors[first_node as usize].insert(second_node);
        backward_neighbors[second_node as usize].insert(first_node);

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

fn decompress_non_terminal_node(x: NodeId, offset: NodeId) -> NodeId {
    if x >= offset {
        (x - offset << 1) + offset
    } else {
        x
    }
}

fn decompress_non_terminals(rules: Rules, offset: NodeId) -> Rules {
    rules
        .into_iter()
        .map(|(k, (v1, v2))| {
            (
                decompress_non_terminal_node(k, offset),
                (
                    decompress_non_terminal_node(v1, offset),
                    decompress_non_terminal_node(v2, offset),
                ),
            )
        })
        .collect()
}

pub fn build_qlines(
    forward_neighbors: &mut Vec<BTreeSet<NodeId>>,
    backward_neighbors: &mut Vec<BTreeSet<NodeId>>,
    digrams: &mut PriorityQueue<(NodeId, NodeId), ColorSet>,
) -> (Rules, NodeId) {
    log::info!("Building qlines for {} digrams", digrams.len());
    println!("D: {:?}", digrams.clone().into_sorted_vec());
    println!("Nf: {:?}", forward_neighbors);
    println!("Nb: {:?}", backward_neighbors);
    let offset = forward_neighbors.len() as NodeId;
    println!("offset: {}", offset);
    let mut rules: Rules = HashMap::new();

    let mut current_max_node_id = offset;

    while digrams.peek().expect("At least one digram").1 .0.len() >= MAX_OCCURENCES {
        let ((u, v), uv_color_set) = digrams.pop().expect("At least one digram");
        println!("Working on digram: {}-{}", u, v);
        println!("\tD: {:?}", digrams.clone().into_sorted_vec());
        println!("\tNf: {:?}", forward_neighbors);
        println!("\tNb: {:?}", backward_neighbors);
        let non_terminal: NodeId = current_max_node_id;
        current_max_node_id += 1;

        backward_neighbors.push(backward_neighbors[u as usize].clone());
        forward_neighbors.push(forward_neighbors[v as usize].clone());

        // println!("{} -> {} | {}", non_terminal, u, v);
        rules.insert(non_terminal, (u, v));
        for w in &backward_neighbors[u as usize] {
            forward_neighbors[*w as usize].insert(non_terminal);

            println!("Working on w: {}", w);
            println!("\tNf: {:?}", forward_neighbors);
            println!("\tNb: {:?}", backward_neighbors);

            let new_wq_set = digrams
                .get_priority(&(*w, u))
                .expect("w-u exists")
                .intersection(&uv_color_set);
            digrams.push((*w, non_terminal), new_wq_set);

            let new_wu_set = digrams
                .get_priority(&(*w, u))
                .expect("w-u exists")
                .difference(&uv_color_set);
            digrams.change_priority(&(*w, u), new_wu_set);
        }
        for w in &forward_neighbors[v as usize] {
            backward_neighbors[*w as usize].insert(non_terminal);

            let new_qw_set = digrams
                .get_priority(&(v, *w))
                .expect("v-w exists")
                .intersection(&uv_color_set);
            digrams.push((non_terminal, *w), new_qw_set);

            let new_vw_set = digrams
                .get_priority(&(v, *w))
                .expect("v-w exists")
                .difference(&uv_color_set);
            digrams.change_priority(&(v, *w), new_vw_set);
        }
        // println!("\t{:?}", digrams.clone().into_sorted_vec());
    }

    // Reintroduce space into non-terminal node ids for reverse non-terminals
    let rules = decompress_non_terminals(rules, offset);

    (rules, offset)
}

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    /// Name of the person to greet
    #[arg(short, long)]
    file: String,
}

fn main() {
    env_logger::init();
    let args = Args::parse();
    let node_ids_by_name = parse_node_ids(&args.file);
    let (mut forward_neighbors, mut backward_neigbors, mut digrams, _path_id_to_path_segment) =
        parse_gfa_paths_walks(&args.file, &node_ids_by_name);
    let (rules, offset) =
        build_qlines(&mut forward_neighbors, &mut backward_neigbors, &mut digrams);
    let mut rules = rules.into_iter().collect::<Vec<_>>();
    rules.sort_by_key(|r| r.0);
    for (target, (x, y)) in rules {
        let x_text = if x >= offset {
            x.to_string()
        } else {
            format!("{}{}", (x >> 1) + 1, if x & 1 == 1 { "+" } else { "-" })
        };
        let y_text = if y >= offset {
            y.to_string()
        } else {
            format!("{}{}", (y >> 1) + 1, if y & 1 == 1 { "+" } else { "-" })
        };
        println!("{} <- {} {}", target, x_text, y_text);
    }
}
