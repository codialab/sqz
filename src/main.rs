mod color_set;
mod compressor;
mod node_id;
mod parser;
mod path_segment;

use clap::Parser;
use color_set::ColorSet;
use compressor::encode_paths;
use node_id::{NodeId, RawNodeId};
use parser::{canonize, parse_gfa_paths_walks, parse_node_ids};
use priority_queue::PriorityQueue;
use std::collections::{HashMap, HashSet};
use std::fmt;

const MAX_OCCURENCES: usize = 2;

type NeighborList = Vec<HashSet<NodeId>>;
type Digrams = PriorityQueue<(NodeId, NodeId), ColorSet>;
type Rules = HashMap<NodeId, Rule>;

#[derive(Clone)]
pub struct Rule {
    left: NodeId,
    right: Vec<NodeId>,
    #[allow(dead_code)]
    colors: ColorSet,
}

impl fmt::Debug for Rule {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{} -> {:?}", self.left, self.right,)
    }
}

impl fmt::Display for Rule {
    // This trait requires `fmt` with this exact signature.
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{} -> {:?}", self.left, self.right,)
    }
}

pub fn build_qlines(neighbors: &mut NeighborList, digrams: &mut Digrams) -> (Rules, NodeId) {
    log::info!("Building qlines for {} digrams", digrams.len());
    let offset = NodeId::from_raw(neighbors.len() as u64);
    let mut rules: Rules = HashMap::new();
    let mut parents: HashMap<NodeId, Option<NodeId>> = HashMap::new();

    let mut current_max_node_id = offset;

    while digrams.peek().expect("At least one digram").1.len() >= MAX_OCCURENCES {
        let ((u, v), uv_color_set) = digrams.pop().expect("At least one digram");
        let non_terminal: NodeId = current_max_node_id;
        current_max_node_id += 2;

        if u >= offset {
            parents
                .entry(u.get_forward())
                .and_modify(|e| *e = None)
                .or_insert(Some(non_terminal));
        }
        if v >= offset {
            parents
                .entry(v.get_forward())
                .and_modify(|e| *e = None)
                .or_insert(Some(non_terminal));
        }

        // Create space in neighbors list
        neighbors.push(HashSet::new());
        neighbors.push(HashSet::new());

        // First store all insertion/deletions that are done later to avoid
        // having to read an mutate neighbors at the same time
        let mut neighbors_to_insert: Vec<(NodeId, NodeId)> = Vec::new();
        let mut neighbors_to_remove: Vec<(NodeId, NodeId)> = Vec::new();

        neighbors_to_remove.push((u, v));
        neighbors_to_remove.push((v.flip(), u.flip()));

        if u == NodeId::new(0, 1) && v == NodeId::new(146, 1) {
            log::error!("==> {} | {} : {:?}, {:?} {:?}", u, v, uv_color_set, neighbors[u.flip().get_idx()], neighbors[v.get_idx()]);
        }

        let mut mutation_outgoing = HashMap::new();
        for n in neighbors.get(u.flip().get_idx()).unwrap() {
            let n = n.flip();
            if n == v {
                continue;
            }

            let nu_set = digrams
                .get_priority(&canonize(n, u))
                .unwrap_or_else(|| {
                    eprintln!(
                        "nu: {} - {} | {:?} | offset: {}",
                        n,
                        u,
                        canonize(n, u),
                        offset
                    );
                    panic!("n-u should exist");
                });

            let is_nu_flipped = canonize(n, u).1 == n;
            let is_nq_flipped = canonize(n, non_terminal).1 == n;
            let (nu_set, nq_set) = uv_color_set.xu_intersection(nu_set, &mut mutation_outgoing, is_nu_flipped, is_nq_flipped);
            if nq_set.is_empty() {
                continue;
            }
            if nu_set.is_empty() {
                neighbors_to_remove.push((n, u));
                neighbors_to_remove.push((u.flip(), n.flip()));
            }

            digrams.push(canonize(n, non_terminal), nq_set);
            digrams.change_priority(&canonize(n, u), nu_set);

            neighbors_to_insert.push((n, non_terminal));
            neighbors_to_insert.push((non_terminal.flip(), n.flip()));
        }
        // for n in &neighbors[v.get_idx()] {
        //     let n = *n;
        //     if n == v && u == v {
        //         // log::debug!("Self loop case: {} -> {} - {} - {}", non_terminal, u, v, n);
        //         neighbors_to_remove.push((v, n));
        //         neighbors_to_remove.push((n.flip(), v.flip()));
        //         continue;
        //     }
        //     if digrams.get_priority(&canonize(v, n)).is_none() {
        //         log::error!("Missing digram: {}-{}, for rule: {}-{}", v, n, u, v);
        //     }
        //     let qn_set = digrams
        //         .get_priority(&canonize(v, n))
        //         .expect("v-n exists")
        //         .intersection(&uv_color_set, 1);
        //     if qn_set.is_empty() {
        //         continue;
        //     }
        //     let vn_set = digrams
        //         .get_priority(&canonize(v, n))
        //         .expect("v-n exists")
        //         .difference(&qn_set);
        //     if canonize(non_terminal, n) == (NodeId::new(0, 1), NodeId::new(146,1)) || canonize(non_terminal, n) == canonize(NodeId::new(137, 0), NodeId::new(0, 0)) {
        //         log::error!("==> q: {} | u: {} | v: {} | n: {} : uv_colorset: {}, orig_vn: {}, qn: {}, new_vn: {}", non_terminal, u, v, n, uv_color_set, digrams.get_priority(&canonize(v, n)).unwrap(), qn_set, vn_set);
        //     }
        //     let mut empty_vn_set = false;
        //     if vn_set.is_empty() {
        //         neighbors_to_remove.push((v, n));
        //         neighbors_to_remove.push((n.flip(), v.flip()));
        //         empty_vn_set = true;
        //     }

        //     if n == u && u != v {
        //         // log::debug!("Loop over rule case: {} -> {} - {} - {}", non_terminal, u, v, n);
        //         digrams.push(canonize(non_terminal, non_terminal), qn_set);
        //         digrams.change_priority(&canonize(v, n), vn_set);

        //         // <n += <(u>v>)
        //         neighbors_to_insert.push((non_terminal.flip(), non_terminal.flip()));
        //         neighbors_to_insert.push((non_terminal, non_terminal));
        //     } else {
        //         digrams.push(canonize(non_terminal, n), qn_set);
        //         digrams.change_priority(&canonize(v, n), vn_set);

        //         // <n += <(u>v>)
        //         neighbors_to_insert.push((n.flip(), non_terminal.flip()));
        //         neighbors_to_insert.push((non_terminal, n));
        //     }

        //     if empty_vn_set {
        //         digrams.remove(&canonize(v, n));
        //     }
        // }
        neighbors_to_insert.into_iter().for_each(|(key, value)| {
            neighbors[key.get_idx()].insert(value);
        });
        neighbors_to_remove.into_iter().for_each(|(key, value)| {
            neighbors[key.get_idx()].remove(&value);
        });

        rules.insert(
            non_terminal,
            Rule {
                left: non_terminal,
                right: vec![u, v],
                colors: uv_color_set,
            },
        );
    }
    log::info!(
        "Built {} rules, more than 0 parents: {}, exactly 1 parent: {}",
        rules.len(),
        parents.len(),
        parents
            .iter()
            .map(|(_, x)| if x.is_some() { 1 } else { 0 })
            .sum::<u64>()
    );

    let rules: Rules = rules
        .into_iter()
        .filter(|(_k, v)| !v.colors.is_empty() || v.colors.len() == 1)
        .collect();
    log::info!("After removal: {}", rules.len());

    let rules = merge_rules(rules, parents, offset);
    log::info!("After merging: {}", rules.len());

    (rules, offset)
}

fn reverse_rule(right: &Vec<NodeId>) -> Vec<NodeId> {
    right.iter().copied().rev().map(|x| x.flip()).collect()
}

fn merge_rules(
    mut rules: Rules,
    parents: HashMap<NodeId, Option<NodeId>>,
    offset: NodeId,
) -> Rules {
    log::info!("Merging rules");
    let mut parents: HashMap<NodeId, NodeId> = parents
        .into_iter()
        .filter_map(|(id, parent)| parent.map(|p| (id, p)))
        .collect();
    let mut mergeables: Vec<NodeId> = parents.keys().copied().collect();
    let mut counter = 0;
    let mut counter_diff_colors = 0;

    log::debug!("Mergeables: {}", mergeables.len());
    while !mergeables.is_empty() {
        let q = mergeables.pop().expect("mergeables should contain a value");
        // log::debug!("Merging: {} (child of {:?}), with rule: {:?} (of rule: {:?})", q, parents.get(&q), rules.get(&q), rules.get(&parents[&q]));
        if rules[&q].colors.equal_set(&rules[&parents[&q]].colors) {
            counter += 1;
            let parent = parents[&q];

            // Insert
            let mut idx = 0;
            let mut reverse = false;
            for (index, node) in rules[&parent].right.iter().enumerate() {
                if node.get_forward() == q {
                    idx = index;
                    reverse = node.get_orientation() == 1;
                }
            }
            let rule_to_insert = if reverse {
                reverse_rule(&rules[&q].right)
            } else {
                rules[&q].right.clone()
            };
            rules
                .get_mut(&parent)
                .expect("rules has parent")
                .right
                .splice(idx..idx + 1, rule_to_insert);

            // Set new parents
            for child in &rules[&q].right {
                if let Some(parent_of_child) = parents.get_mut(&child.get_forward()) {
                    *parent_of_child = parent;
                }
            }

            rules.remove(&q);
        } else {
            counter_diff_colors += 1;
        }
    }
    log::info!(
        "Merged {} rules, {} times the color set was different",
        counter,
        counter_diff_colors,
    );
    rules
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
    let (rules, offset) = build_qlines(&mut neighbors, &mut digrams);
    let encoded_paths = encode_paths(&args.file, &rules, offset, &node_ids_by_name);
    let rules = rules.into_iter().collect::<Vec<_>>();
    for (_, rule) in rules {
        println!("Q\t{} | {:?}", rule, rule.colors);
    }
    for (path_name, path) in encoded_paths {
        println!("Z\t{}\t{:?}", path_name, path);
    }
}
