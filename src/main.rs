mod path_segment;
mod parser;

use clap::Parser;
use priority_queue::PriorityQueue;
use std::collections::{HashMap, HashSet};
use parser::{parse_node_ids, parse_gfa_paths_walks, canonize};
use std::fmt;
use std::ops::{Add, AddAssign, Sub};

const MAX_OCCURENCES: usize = 2;

type RawNodeId = u64;

#[derive(Clone, Copy, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub struct NodeId(u64);

impl NodeId {
    fn new(id: RawNodeId, orientation: u64) -> Self {
        Self(id << 1 ^ orientation)
    }
    fn flip(&self) -> NodeId {
        Self(self.0 ^ 1)
    }

    fn get_id(&self) -> RawNodeId {
        self.0 >> 1
    }

    fn get_orientation(&self) -> u64 {
        self.0 & 1
    }

    fn get_idx(&self) -> usize {
        self.0 as usize
    }
}

impl Add<u64> for NodeId {
    type Output = Self;

    fn add(self, other: u64) -> Self {
        Self(self.0 + other)
    }
}

impl AddAssign<u64> for NodeId {
    fn add_assign(&mut self, other: u64) {
        *self = Self(self.0 + other);
    }
}

impl Sub for NodeId {
    type Output = u64;

    fn sub(self, other: Self) -> Self::Output {
        self.0 - other.0
    }
}

impl fmt::Debug for NodeId {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}{}", if self.0 & 1 == 0 { '>' } else { '<' }, (self.0 >> 1) + 1)
    }
}

impl fmt::Display for NodeId {
    // This trait requires `fmt` with this exact signature.
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}{}", if self.0 & 1 == 0 { '>' } else { '<' }, (self.0 >> 1) + 1)
    }
}

type NeighborList = Vec<HashSet<NodeId>>;
type Digrams = PriorityQueue<(NodeId, NodeId), ColorSet>;
type Rules = HashMap<NodeId, (NodeId, NodeId)>;


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

    fn is_empty(&self) -> bool {
        self.0.is_empty()
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


fn decompress_non_terminal_node(x: NodeId, offset: NodeId) -> NodeId {
    if x >= offset {
        offset + ((x - offset) << 1)
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
    neighbors: &mut NeighborList,
    digrams: &mut Digrams,
) -> (Rules, NodeId) {
    log::info!("Building qlines for {} digrams", digrams.len());
    println!("D: {:?}", digrams.clone().into_sorted_vec());
    println!("N: {:?}", neighbors);
    let offset = NodeId(neighbors.len() as u64);
    println!("offset: {}", offset);
    let mut rules: Rules = HashMap::new();

    let mut current_max_node_id = offset;

    while digrams.peek().expect("At least one digram").1 .0.len() >= MAX_OCCURENCES {
        let ((u, v), uv_color_set) = digrams.pop().expect("At least one digram");
        let non_terminal: NodeId = current_max_node_id;
        println!("Working on digram: {}, {} <- {} | {}", non_terminal, non_terminal.0, u, v);
        println!("\tD: {:?}", digrams.clone().into_sorted_vec());
        println!("\tN: {:?}", neighbors);
        current_max_node_id += 2;

        rules.insert(non_terminal, (u, v));

        // Create space in neighbors list
        neighbors.push(HashSet::new());
        neighbors.push(HashSet::new());

        // First store all insertion/deletions that are done later to avoid
        // having to read an mutate neighbors at the same time
        let mut neighbors_to_insert: Vec<(NodeId, NodeId)> = Vec::new();
        let mut neighbors_to_remove: Vec<(NodeId, NodeId)> = Vec::new();
        for n in &neighbors[v.get_idx()] {
            let n = *n;
            println!("vn: {} - {} | {:?}", v, n, canonize(v, n));
            let qn_set = digrams.get_priority(&canonize(v, n)).expect("v-n exists").intersection(&uv_color_set);
            if qn_set.is_empty() {
                continue;
            }
            digrams.push(canonize(non_terminal, n), qn_set);

            let vn_set = digrams.get_priority(&canonize(v, n)).expect("v-n exists").difference(&uv_color_set);
            if vn_set.is_empty() {
                neighbors_to_remove.push((v, n));
                neighbors_to_remove.push((n.flip(), v.flip()));
            }
            digrams.change_priority(&canonize(v, n), vn_set);

            // <n += <(u>v>)
            neighbors_to_insert.push((n.flip(), non_terminal.flip()));
            neighbors_to_insert.push((non_terminal, n));
        }
        for n in neighbors.get(u.flip().get_idx()).unwrap() {
            let n = n.flip();
            println!("nu: {} - {} | {:?}", n, u, canonize(n, u));
            let nq_set = digrams.get_priority(&canonize(n, u)).expect("n-u exists").intersection(&uv_color_set);
            if nq_set.is_empty() {
                continue;
            }
            digrams.push(canonize(n, non_terminal), nq_set);

            let nu_set = digrams.get_priority(&canonize(n, u)).expect("n-u exists").difference(&uv_color_set);
            if nu_set.is_empty() {
                neighbors_to_remove.push((n, u));
                neighbors_to_remove.push((u.flip(), n.flip()));
            }
            digrams.change_priority(&canonize(n, u), nu_set);

            neighbors_to_insert.push((n, non_terminal));
            neighbors_to_insert.push((non_terminal.flip(), n.flip()));
            println!("\tnti: {:?}", neighbors_to_insert);
        }
        neighbors_to_insert.into_iter().for_each(|(key, value)| {
            neighbors[key.get_idx()].insert(value);
        });
        neighbors_to_remove.into_iter().for_each(|(key, value)| {
            //let idx = neighbors[key as usize].iter().position(|x| *x == value).expect("can remove neighbor");
            neighbors[key.get_idx()].remove(&value);
        });
        // println!("\t{:?}", digrams.clone().into_sorted_vec());

        println!("\tFinished D: {:?}", digrams.clone().into_sorted_vec());
    }


    println!("Final N: {:?}", neighbors);
    // Reintroduce space into non-terminal node ids for reverse non-terminals
    //let rules = decompress_non_terminals(rules, offset);

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
    let (mut neighbors, mut digrams, _path_id_to_path_segment) =
        parse_gfa_paths_walks(&args.file, &node_ids_by_name);
    let (rules, _offset) =
        build_qlines(&mut neighbors, &mut digrams);
    let mut rules = rules.into_iter().collect::<Vec<_>>();
    rules.sort_by_key(|r| r.0);
    for (target, (x, y)) in rules {
        println!("{} <- {} {}", target, x, y);
    }
}
