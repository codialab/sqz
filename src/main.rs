mod path_segment;
mod parser;

use clap::Parser;
use priority_queue::PriorityQueue;
use std::collections::{BTreeSet, HashMap, HashSet};
use parser::{parse_node_ids, parse_gfa_paths_walks, canonize};

const MAX_OCCURENCES: usize = 2;

type NodeId = u64;
type NeighborList = Vec<BTreeSet<NodeId>>;
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

    fn is_disjoint(&self, other: &ColorSet) -> bool {
        self.0.is_disjoint(&other.0)
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
        ((x - offset) << 1) + offset
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

fn flip(x: NodeId) -> NodeId {
    x ^ 1
}

pub fn build_qlines(
    neighbors: &mut Vec<BTreeSet<NodeId>>,
    digrams: &mut PriorityQueue<(NodeId, NodeId), ColorSet>,
) -> (Rules, NodeId) {
    log::info!("Building qlines for {} digrams", digrams.len());
    println!("D: {:?}", digrams.clone().into_sorted_vec());
    println!("N: {:?}", neighbors);
    let offset = neighbors.len() as NodeId;
    println!("offset: {}", offset);
    let mut rules: Rules = HashMap::new();

    let mut current_max_node_id = offset;

    while digrams.peek().expect("At least one digram").1 .0.len() >= MAX_OCCURENCES {
        let ((u, v), uv_color_set) = digrams.pop().expect("At least one digram");
        println!("Working on digram: {}-{}", u, v);
        println!("\tD: {:?}", digrams.clone().into_sorted_vec());
        println!("\tN: {:?}", neighbors);
        let non_terminal: NodeId = current_max_node_id;
        current_max_node_id += 2;

        // neighbors.push(neighbors[v as usize].iter()
        //                .filter(|n| !digrams.get_priority(&canonize(v, n)).expect("v-n exists").is_disjoint(&uv_color_set)).collect::<BTreeSet<_>>());
        // neighbors.push(neighbors[flip(u) as usize].clone());

        // println!("{} -> {} | {}", non_terminal, u, v);
        rules.insert(non_terminal, (u, v));
        for n in &neighbors[v as usize] {
            // TODO check if invalid neighbor (i.e. intersection is empty)

            // <n += <(u>v>)
            neighbors[flip(*n) as usize].insert(flip(non_terminal));

            println!("Working on w: {}", w);
            println!("\tNf: {:?}", neighbors);

            let new_qw_set = digrams
                .get_priority(&canonize(*w, u))
                .expect("w-u exists")
                .intersection(&uv_color_set);
            digrams.push((*w, non_terminal), new_wq_set);

            let new_wu_set = digrams
                .get_priority(&(*w, u))
                .expect("w-u exists")
                .difference(&uv_color_set);
            digrams.change_priority(&(*w, u), new_wu_set);
        }
        for w in &neighbors[flip(u) as usize] {
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
    let (mut neighbors, mut digrams, _path_id_to_path_segment) =
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
