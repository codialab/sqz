mod color_set;
mod node_id;
mod parser;
mod path_segment;

use clap::Parser;
use color_set::ColorSet;
use node_id::{NodeId, RawNodeId};
use parser::{canonize, parse_gfa_paths_walks, parse_node_ids};
use priority_queue::PriorityQueue;
use std::collections::{HashMap, HashSet};
use std::fmt;

const MAX_OCCURENCES: usize = 2;

type NeighborList = Vec<HashSet<NodeId>>;
type Digrams = PriorityQueue<(NodeId, NodeId), ColorSet>;
type Rules = HashMap<NodeId, (NodeId, NodeId)>;

struct Rule {
    left: NodeId,
    right: Vec<NodeId>,
    colors: ColorSet,
}

impl fmt::Debug for Rule {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "{} -> {:?}", self.left, self.right,
        )
    }
}

impl fmt::Display for Rule {
    // This trait requires `fmt` with this exact signature.
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "{} -> {:?}", self.left, self.right,
        )
    }
}

pub fn build_qlines(neighbors: &mut NeighborList, digrams: &mut Digrams) -> (Rules, NodeId) {
    log::info!("Building qlines for {} digrams", digrams.len());
    // println!("D: {:?}", digrams.clone().into_sorted_vec());
    // println!("N: {:?}", neighbors);
    let offset = NodeId::from_raw(neighbors.len() as u64);
    // println!("offset: {}", offset);
    let mut rules: Rules = HashMap::new();

    let mut current_max_node_id = offset;

    while digrams.peek().expect("At least one digram").1.len() >= MAX_OCCURENCES {
        let ((u, v), uv_color_set) = digrams.pop().expect("At least one digram");
        let non_terminal: NodeId = current_max_node_id;
        println!(
            "Working on digram: {}, {} <- {} | {}",
            non_terminal,
            non_terminal.get_idx(),
            u,
            v
        );
        // println!("\tD: {:?}", digrams.clone().into_sorted_vec());
        // println!("\tN: {:?}", neighbors);
        current_max_node_id += 2;

        rules.insert(non_terminal, (u, v));

        // Create space in neighbors list
        neighbors.push(HashSet::new());
        neighbors.push(HashSet::new());

        // First store all insertion/deletions that are done later to avoid
        // having to read an mutate neighbors at the same time
        let mut neighbors_to_insert: Vec<(NodeId, NodeId)> = Vec::new();
        let mut neighbors_to_remove: Vec<(NodeId, NodeId)> = Vec::new();

        neighbors_to_remove.push((u, v));
        neighbors_to_remove.push((v.flip(), u.flip()));

        for n in &neighbors[v.get_idx()] {
            let n = *n;
            // println!("vn: {} - {} | {:?}", v, n, canonize(v, n));
            let qn_set = digrams
                .get_priority(&canonize(v, n))
                .expect("v-n exists")
                .intersection(&uv_color_set, 1);
            if qn_set.is_empty() {
                continue;
            }
            let vn_set = digrams
                .get_priority(&canonize(v, n))
                .expect("v-n exists")
                .difference(&qn_set);
            if vn_set.is_empty() {
                neighbors_to_remove.push((v, n));
                neighbors_to_remove.push((n.flip(), v.flip()));
            }

            digrams.push(canonize(non_terminal, n), qn_set);
            digrams.change_priority(&canonize(v, n), vn_set);

            // <n += <(u>v>)
            neighbors_to_insert.push((n.flip(), non_terminal.flip()));
            neighbors_to_insert.push((non_terminal, n));
        }
        for n in neighbors.get(u.flip().get_idx()).unwrap() {
            let n = n.flip();
            let nq_set = digrams
                .get_priority(&canonize(n, u))
                .unwrap_or_else(|| {
                    println!("nu: {} - {} | {:?} | offset: {}", n, u, canonize(n, u), offset);
                    panic!("n-u should exist");
                })
                .intersection(&uv_color_set, 0);
            if nq_set.is_empty() {
                continue;
            }
            let nu_set = digrams
                .get_priority(&canonize(n, u))
                .expect("n-u exists")
                .difference(&nq_set);
            if nu_set.is_empty() {
                neighbors_to_remove.push((n, u));
                neighbors_to_remove.push((u.flip(), n.flip()));
            }

            digrams.push(canonize(n, non_terminal), nq_set);
            digrams.change_priority(&canonize(n, u), nu_set);

            neighbors_to_insert.push((n, non_terminal));
            neighbors_to_insert.push((non_terminal.flip(), n.flip()));
            // println!("\tnti: {:?}", neighbors_to_insert);
        }
        neighbors_to_insert.into_iter().for_each(|(key, value)| {
            neighbors[key.get_idx()].insert(value);
        });
        neighbors_to_remove.into_iter().for_each(|(key, value)| {
            //let idx = neighbors[key as usize].iter().position(|x| *x == value).expect("can remove neighbor");
            neighbors[key.get_idx()].remove(&value);
        });
        // println!("\t{:?}", digrams.clone().into_sorted_vec());

        // println!("\tFinished D: {:?}", digrams.clone().into_sorted_vec());
    }

    // println!("Final N: {:?}", neighbors);
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
    let (rules, _offset) = build_qlines(&mut neighbors, &mut digrams);
    let mut rules = rules.into_iter().collect::<Vec<_>>();
    rules.sort_by_key(|r| r.0);
    for (target, (x, y)) in rules {
        println!("{} <- {} {}", target, x, y);
    }
}
