mod color_set;
mod compressor;
mod node_id;
mod parser;
mod path_segment;

use clap::Parser;
use color_set::OrdColorSet;
use compressor::encode_paths2;
use indexmap::IndexMap;
use node_id::{NodeId, RawNodeId};
use parser::{canonize, parse_gfa_paths_walks, parse_node_ids};
use priority_queue::PriorityQueue;
use std::collections::{HashMap, HashSet};
use std::fmt;
use std::mem;

const MAX_OCCURENCES: usize = 2;

type NeighborList = Vec<HashSet<NodeId>>;
type Digrams = PriorityQueue<(NodeId, NodeId), OrdColorSet>;
type Rules = IndexMap<NodeId, Rule>;

#[derive(Clone)]
pub struct Rule {
    left: NodeId,
    right: Vec<NodeId>,
    colors: OrdColorSet,
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
    let mut rules: Rules = IndexMap::new();
    let mut parents: IndexMap<NodeId, Option<NodeId>> = IndexMap::new();

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

        insert_edge(&mut neighbors_to_remove, u, v);

        let mut mutation_outgoing = HashMap::new();

        let mut new_uv_set = None;

        // let should_print = non_terminal == NodeId::new(11, 0) || non_terminal == NodeId::new(8, 0) || non_terminal == NodeId::new(9, 0) || non_terminal == NodeId::new(10, 0); //(u.get_forward() == NodeId::new(9, 0) && v.get_forward() == NodeId::new(991, 0)) || (u.get_forward() == NodeId::new(987, 0) && v.get_forward() == NodeId::new(989, 0));
        let should_print = true;

        for n in neighbors.get(u.flip().get_idx()).unwrap() {
            let n = n.flip();
            println!("nu's n: {}", n);

            if n == u && u == v {
                insert_edge(&mut neighbors_to_remove, n, u);
                continue;
            }

            let nu_set = digrams.get_priority(&canonize(n, u)).unwrap_or_else(|| {
                log::error!(
                    "nu: {} - {} | {:?} | offset: {}",
                    n,
                    u,
                    canonize(n, u),
                    offset
                );
                panic!("n-u should exist");
            });

            let is_nu_flipped = is_edge_flipped(n, u);
            let is_nq_flipped = is_edge_flipped(n, non_terminal);
            let (new_nu_set, mut nq_set) = uv_color_set.colors.xu_intersection(
                &nu_set.colors,
                &mut mutation_outgoing,
                is_nu_flipped,
                is_nq_flipped,
            );
            if should_print {
                println!("=====================");
                println!("n: {}, u: {}, v: {}, (q: {})", n, u, v, non_terminal);
                println!(
                    "nu: {:?}, uv: {:?}, nu_flipped: {}, nq_flipped: {}",
                    nu_set, uv_color_set, is_nu_flipped, is_nq_flipped
                );
                println!(
                    "nu: {:?}, nq: {:?}, mutation: {:?}",
                    new_nu_set, nq_set, mutation_outgoing
                );
                println!("=====================");
            }
            if nq_set.is_empty() {
                continue;
            }
            if new_nu_set.is_empty() {
                insert_edge(&mut neighbors_to_remove, n, u);
            }

            if n == v && u != v {
                if should_print {
                    println!("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
                    println!("non_terminal: {}, nq_set: {:?}, mut: {:?}", non_terminal, nq_set, mutation_outgoing);
                    println!("uv: {:?}, new_nq: {:?}", uv_color_set, nq_set);
                    println!("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
                }
                let real_nq_set = nq_set.cleanup_pre_self_loop(&uv_color_set.colors, is_nq_flipped);
                if is_nq_flipped == is_edge_flipped(non_terminal, non_terminal) {
                    digrams.push(
                        canonize(non_terminal, non_terminal),
                        OrdColorSet::new(nq_set, true),
                    );
                } else {
                    nq_set.flip_all();
                    digrams.push(
                        canonize(non_terminal, non_terminal),
                        OrdColorSet::new(nq_set, true),
                    );
                }
                if !real_nq_set.is_empty() {
                    digrams.push(
                        canonize(n, non_terminal),
                        OrdColorSet::new(real_nq_set, false),
                    );
                    insert_edge(&mut neighbors_to_insert, n, non_terminal);
                }
                //new_nu_set.add_addition(nu_set_addition);
                digrams.change_priority(&canonize(n, u), OrdColorSet::new(new_nu_set, false));

                insert_edge(&mut neighbors_to_insert, non_terminal, non_terminal);
            } else {
                digrams.push(canonize(n, non_terminal), OrdColorSet::new(nq_set, false));
                digrams.change_priority(&canonize(n, u), OrdColorSet::new(new_nu_set, n == u));

                insert_edge(&mut neighbors_to_insert, n, non_terminal);
            }
        }

        let mut self_sets = if u == v {
            log::debug!("Self loop case: {} -> {} - {}", non_terminal, u, v);
            let self_sets = uv_color_set.colors.sectionize(&mut mutation_outgoing);
            // log::error!("Self loop contents: qq: {:?}, qv: {:?}, uv: {:?}", self_sets.0, self_sets.1, self_sets.2);
            Some(self_sets)
        } else {
            None
        };

        for n in &neighbors[v.get_idx()] {
            let n = *n;
            println!("vn's n: {}", n);

            if n == u && u != v {
                // TODO: handle this case
                let vn_set = digrams.get_priority(&canonize(v, n)).unwrap_or_else(|| {
                    log::error!(
                        "vn: {} - {} | {:?} | offset: {}",
                        v,
                        n,
                        canonize(v, n),
                        offset
                    );
                    panic!("v-n should exist");
                });
                // let qq_set = digrams.get_priority(&canonize(non_terminal, non_terminal)).unwrap_or_else(|| {
                //     log::error!(
                //         "qq: {} - {} | {:?} | offset: {}",
                //         non_terminal,
                //         non_terminal,
                //         canonize(non_terminal, non_terminal),
                //         offset
                //     );
                //     panic!("q-q should exist");
                // });
                println!("++++++++++++++++++++++++++");
                println!("Pre-self-loop case 2: {} -> {} {} (n: {})", non_terminal, u, v, n);
                println!("vn: {:?}", vn_set);
                // println!("qq: {:?}", qq_set);
                println!("uv: {:?}", uv_color_set);
                println!("++++++++++++++++++++++++++");
            } else if n == v && u == v {
                let (mut qq_set, mut qv_set, uv_temp) =
                    mem::take(&mut self_sets).expect("self sets should have been set");
                println!("qq_set: {:?}", qq_set);
                new_uv_set = Some(uv_temp);
                if !qq_set.is_empty() {
                    if is_edge_flipped(non_terminal, non_terminal) {
                        qq_set.flip_all();
                        digrams.push(
                            canonize(non_terminal, non_terminal),
                            OrdColorSet::new(qq_set, true),
                        );
                    } else {
                        digrams.push(
                            canonize(non_terminal, non_terminal),
                            OrdColorSet::new(qq_set, true),
                        );
                    }
                    insert_edge(&mut neighbors_to_insert, non_terminal, non_terminal);
                }
                if !qv_set.is_empty() {
                    if is_edge_flipped(non_terminal, v) {
                        qv_set.flip_all();
                        digrams.push(canonize(non_terminal, v), OrdColorSet::new(qv_set, false));
                    } else {
                        digrams.push(canonize(non_terminal, v), OrdColorSet::new(qv_set, false));
                    }
                    insert_edge(&mut neighbors_to_insert, non_terminal, v);
                }
                insert_edge(&mut neighbors_to_remove, v, n);
                continue;
            }

            let vn_set = digrams.get_priority(&canonize(v, n)).unwrap_or_else(|| {
                log::error!(
                    "vn: {} - {} | {:?} | offset: {}",
                    v,
                    n,
                    canonize(v, n),
                    offset
                );
                panic!("v-n should exist");
            });
            let is_vn_flipped = is_edge_flipped(v, n);
            let is_qn_flipped = is_edge_flipped(non_terminal, n);
            let (mut new_vn_set, mut qn_set) = if u != v {
                uv_color_set.colors.vy_intersection(
                    &vn_set.colors,
                    true,
                    is_vn_flipped,
                    is_qn_flipped,
                )
            } else {
                // let dont_mutate = HashMap::new();
                uv_color_set.colors.vy_intersection(
                    &vn_set.colors,
                    false,
                    is_vn_flipped,
                    is_qn_flipped,
                )
            };
            if should_print {
                println!("=====================");
                println!("u: {}, v: {}, (q: {}), n: {}", u, v, non_terminal, n);
                println!(
                    "vn: {:?}, uv: {:?}, vn_flipped: {}, qn_flipped: {}",
                    vn_set, uv_color_set, is_vn_flipped, is_qn_flipped
                );
                println!(
                    "vn: {:?}, qn: {:?}, mutation: {:?}",
                    new_vn_set, qn_set, mutation_outgoing
                );
                println!("=====================");
            }

            // Reduce qn_set further to only include the edge only in the case that it was an odd number of self-loops
            if u == v {
                log::debug!("Running on qn_set: {} -> {} {}, n: {}, old_qn: {:?}, mutation: {:?}", non_terminal, u, v, n, qn_set, mutation_outgoing);
                let (new_qn_set, vn_set_addition) =
                    qn_set.self_vy_intersection(&mutation_outgoing, is_vn_flipped, is_qn_flipped);
                qn_set = new_qn_set;
                new_vn_set.add_addition(vn_set_addition);
                log::debug!("Running on qn_set: qn: {:?}, new_vn_set: {:?}", qn_set, new_vn_set);
            }

            if qn_set.is_empty() {
                continue;
            }

            let mut empty_vn_set = false;
            if new_vn_set.is_empty() {
                insert_edge(&mut neighbors_to_remove, v, n);
                empty_vn_set = true;
            }

            if is_qn_flipped == is_vn_flipped {
                digrams.push(canonize(non_terminal, n), OrdColorSet::new(qn_set, false));
            } else {
                // qn_set.flip_all();
                digrams.push(canonize(non_terminal, n), OrdColorSet::new(qn_set, false));
            }
            digrams.change_priority(&canonize(v, n), OrdColorSet::new(new_vn_set, v == n));

            insert_edge(&mut neighbors_to_insert, non_terminal, n);

            if empty_vn_set {
                digrams.remove(&canonize(v, n));
            }
        }
        neighbors_to_insert.into_iter().for_each(|(key, value)| {
            neighbors[key.get_idx()].insert(value);
        });
        neighbors_to_remove.into_iter().for_each(|(key, value)| {
            neighbors[key.get_idx()].remove(&value);
        });

        if new_uv_set.is_none() {
            rules.insert(
                non_terminal,
                Rule {
                    left: non_terminal,
                    right: vec![u, v],
                    colors: uv_color_set,
                },
            );
        } else {
            println!("Self sets: {:?}", self_sets);
            rules.insert(
                non_terminal,
                Rule {
                    left: non_terminal,
                    right: vec![u, v],
                    colors: OrdColorSet::new(new_uv_set.unwrap(), true),
                },
            );
        }
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
        .filter(|(_k, v)| !v.colors.colors.is_empty() || v.colors.len() == 1)
        .collect();
    log::info!("After removal: {}", rules.len());

    // let rules = merge_rules(rules, parents);
    // log::info!("After merging: {}", rules.len());

    (rules, offset)
}

fn is_edge_flipped(a: NodeId, b: NodeId) -> bool {
    if a == NodeId::new(74, 0) && b == NodeId::new(26, 0) {
        println!(
            "FLIP?: {:?}, {:?}, {}, {}",
            (a, b),
            canonize(a, b),
            canonize(a, b).0,
            canonize(a, b).0 == b.flip()
        )
    }
    canonize(a, b).0 == b.flip()
}

fn insert_edge(neighbors_to_insert: &mut Vec<(NodeId, NodeId)>, a: NodeId, b: NodeId) {
    neighbors_to_insert.push((a, b));
    neighbors_to_insert.push((b.flip(), a.flip()));
}

fn reverse_rule(right: &Vec<NodeId>) -> Vec<NodeId> {
    right.iter().copied().rev().map(|x| x.flip()).collect()
}

fn merge_rules(mut rules: Rules, parents: IndexMap<NodeId, Option<NodeId>>) -> Rules {
    log::info!("Merging rules");
    let mut parents: IndexMap<NodeId, NodeId> = parents
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
        if rules[&q]
            .colors
            .colors
            .equal_set(&rules[&parents[&q]].colors.colors)
        {
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
    for (k, v) in digrams {
        println!("digram: {:?}, {:?}", k, v.colors);
    }
    let encoded_paths = encode_paths2(&args.file, &rules, offset, &node_ids_by_name);
    let mut rules = rules.into_iter().collect::<Vec<_>>();
    rules.sort_by_key(|r| r.0.get_idx());
    for (_, rule) in rules {
        println!("Q\t{} | {:?}", rule, rule.colors);
    }
    for (path_name, path) in encoded_paths {
        println!("Z\t{}\t{:?}", path_name, path);
    }
}
