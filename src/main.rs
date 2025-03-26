mod path_segment;
mod parser;

use clap::Parser;
use priority_queue::PriorityQueue;
use std::collections::{HashMap, HashSet};
use parser::{parse_node_ids, parse_gfa_paths_walks, canonize};

const MAX_OCCURENCES: usize = 2;

type NodeId = u64;
type NeighborList = Vec<Vec<NodeId>>;
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
    neighbors: &mut NeighborList,
    digrams: &mut Digrams,
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

        rules.insert(non_terminal, (u, v));

        // First store all insertion/deletions that are done later to avoid
        // having to read an mutate neighbors at the same time
        let mut neighbors_to_insert: Vec<(NodeId, NodeId)> = Vec::new();
        let mut neighbors_to_remove: Vec<(NodeId, NodeId)> = Vec::new();
        for n in &neighbors[v as usize] {
            let qn_set = digrams.get_priority(&canonize(v, *n)).expect("v-n exists").intersection(&uv_color_set);
            if qn_set.is_empty() {
                continue;
            }
            digrams.push(canonize(non_terminal, *n), qn_set);

            let vn_set = digrams.get_priority(&canonize(v, *n)).expect("v-n exists").difference(&uv_color_set);
            if vn_set.is_empty() {
                neighbors_to_remove.push((v, *n));
                neighbors_to_remove.push((flip(*n), flip(v)));
            }
            digrams.change_priority(&canonize(v, *n), vn_set);

            // <n += <(u>v>)
            neighbors_to_insert.push((flip(*n), flip(non_terminal)));
            neighbors_to_insert.push((non_terminal, *n));
        }
        for n in neighbors.get(flip(u) as usize).unwrap() {
            let nq_set = digrams.get_priority(&canonize(*n, u)).expect("n-u exists").intersection(&uv_color_set);
            if nq_set.is_empty() {
                continue;
            }
            digrams.push(canonize(*n, non_terminal), nq_set);

            let nu_set = digrams.get_priority(&canonize(*n, u)).expect("n-u exists").difference(&uv_color_set);
            if nu_set.is_empty() {
                neighbors_to_remove.push((*n, u));
                neighbors_to_remove.push((flip(u), flip(*n)));
            }
            digrams.change_priority(&canonize(*n, u), nu_set);

            neighbors_to_insert.push((*n, non_terminal));
            neighbors_to_insert.push((flip(non_terminal), flip(*n)));
        }
        neighbors_to_insert.into_iter().for_each(|(key, value)| {
            neighbors[key as usize].push(value);
        });
        neighbors_to_remove.into_iter().for_each(|(key, value)| {
            let idx = neighbors[key as usize].iter().position(|x| *x == value).expect("can remove neighbor");
            neighbors[key as usize].swap_remove(idx);
        });
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
        build_qlines(&mut neighbors, &mut digrams);
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
