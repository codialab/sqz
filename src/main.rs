mod path_segment;

use core::str;
use std::{collections::{BinaryHeap, HashMap, HashSet, BTreeSet, BTreeMap}, io::{BufRead, BufReader, Read}};
use clap::Parser;
use flate2::read::MultiGzDecoder;
use itertools::Itertools;
use path_segment::PathSegment;
use std::str::FromStr;

type NodeId = u64;
type PathId = u64;

#[derive(PartialEq, Eq)]
pub struct ColorSet(HashSet<u64>);

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

pub fn parse_gfa_paths_walks(
    gfa_file: &str,
    node_ids_by_name: &HashMap<Vec<u8>, NodeId>,
) -> () {
    log::info!("parsing path + walk sequences");
    let mut data = bufreader_from_compressed_gfa(gfa_file);

    let mut forward_neighbors: Vec<Vec<u64>> = Vec::new();
    let mut backward_neighbors: Vec<Vec<u64>> = Vec::new();
    let mut digrams: BinaryHeap<ColorSet> = BinaryHeap::new();

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
                b'P' => parse_path_seq(&buf, path_id, node_ids_by_name, &mut forward_neighbors, &mut backward_neighbors, &mut digrams),
                b'W' => parse_walk_seq(&buf, node_ids_by_name, &mut forward_neighbors, &mut backward_neighbors, &mut digrams),
                _ => unreachable!(),
            }
            path_id += 1;
        }
    }
}

pub fn parse_path_seq(
    data: &[u8],
    path_id: u64,
    node_ids_by_name: &HashMap<Vec<u8>, NodeId>,
    forward_neighbors: &mut HashMap<NodeId, BTreeSet<NodeId>>,
    backward_neighbors: &mut HashMap<NodeId, BTreeSet<NodeId>>,
    digrams: &mut BTreeMap<(NodeId, NodeId), ColorSet>,
) {
    let mut nodes_visited: HashMap<NodeId, usize> = HashMap::new();
    let mut it = data.iter();
    let end = it
        .position(|x| x == &b'\t' || x == &b'\n' || x == &b'\r')
        .unwrap();

    let prev_node = data[..end]
        .split(|&x| x == b',')
        .take(1)
        .collect::<Vec<_>>()[0]
        .trim_ascii();
    let mut prev_node = node_ids_by_name[&prev_node[..prev_node.len()-1]] | (((prev_node[prev_node.len()-1] == b'+') as NodeId) << (NodeId::BITS - 1));
    data[..end].split(|&x| x == b',').for_each(|current_node| {

        let current_node = current_node.trim_ascii();
        let current_node = node_ids_by_name[&current_node[..current_node.len()-1]] | (((current_node[current_node.len()-1] == b'+') as NodeId) << (NodeId::BITS - 1));

        nodes_visited.entry(current_node).and_modify(|counter| *counter += 1).or_insert(0);
        let count = nodes_visited[&current_node];
        let current_node = (count as u64) << 32 ^ current_node;

        forward_neighbors.entry(prev_node).and_modify(|n| { n.insert(current_node); } ).or_insert(BTreeSet::from([current_node]));
        backward_neighbors.entry(current_node).and_modify(|n| { n.insert(prev_node); }).or_insert(BTreeSet::from([prev_node]));

    });

    log::debug!("parsing path sequences of size {} bytes..", end);
}


pub fn parse_walk_seq(
    data: &[u8],
    node_ids_by_name: &HashMap<Vec<u8>, NodeId>,
    forward_neighbors: &mut Vec<Vec<u64>>,
    backward_neighbors: &mut Vec<Vec<u64>>,
    digrams: &mut BinaryHeap<ColorSet>,
) {
    let mut nodes_visited: HashMap<PathSegment, usize> = HashMap::new();
    let mut it = data.iter();
    let end = it
        .position(|x| x == &b'\t' || x == &b'\n' || x == &b'\r')
        .unwrap();

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

pub fn parse_node_ids(
    gfa_file: &str,
) -> HashMap<Vec<u8>, u64> {
    let mut node2id: HashMap<Vec<u8>, NodeId> = HashMap::default();

    log::info!("constructing indexes for node/edge IDs, node lengths, and P/W lines..");
    let mut node_id = 1; // important: id must be > 0, otherwise counting procedure will produce errors

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
        }         buf.clear();
    }

    log::info!(
        "found: {} nodes",
        node2id.len()
    );
    node2id
}

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    /// Name of the person to greet
    #[arg(short, long)]
    file: String,
}

fn main() {
    let args = Args::parse();
    let node_ids_by_name = parse_node_ids(&args.file);
    parse_gfa_paths_walks(&args.file, &node_ids_by_name);
    println!("Hello, world!");
}
