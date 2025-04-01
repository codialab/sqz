mod color_set;
mod node_id;
mod parser;
mod path_segment;

use clap::Parser;
use color_set::ColorSet;
use node_id::{NodeId, RawNodeId};
use parser::{canonize, parse_gfa_paths_walks, get_nodes_path, get_nodes_walk, parse_node_ids, parse_path_identifier, parse_walk_identifier, bufreader_from_compressed_gfa};
use priority_queue::PriorityQueue;
use std::collections::{HashMap, HashSet, VecDeque};
use std::fmt;
use std::io::BufRead;

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

pub fn encode_paths(file: &str, rules: &Rules, offset: NodeId, node_ids_by_name: &HashMap<Vec<u8>, RawNodeId>) -> HashMap<String, Vec<NodeId>> {
    let mut result: HashMap<String, Vec<NodeId>> = HashMap::new();
    let rules: HashMap<Vec<NodeId>, Rule> = rules.iter().map(|(_k, v)| (v.right.clone(), v.clone())).collect();
    let mut data = bufreader_from_compressed_gfa(file);
    let mut statistics: HashMap<Vec<NodeId>, ColorSet> = rules.iter().map(|(k, _v)| (k.clone(), ColorSet::new())).collect();
    let mut statistics2: HashMap<Vec<NodeId>, usize> = rules.iter().map(|(k, _v)| (k.clone(), 0)).collect();
    let mut path_id = 0;

    let mut buf = vec![];
    while data.read_until(b'\n', &mut buf).unwrap_or(0) > 0 {
        if buf[0] == b'P' || buf[0] == b'W' {
            let (path_seg, buf_path_seg) = match buf[0] {
                b'P' => parse_path_identifier(&buf),
                b'W' => parse_walk_identifier(&buf),
                _ => unreachable!(),
            };

            let nodes = match buf[0] {
                b'P' => get_nodes_path(buf_path_seg, node_ids_by_name),
                b'W' => get_nodes_walk(buf_path_seg, node_ids_by_name),
                _ => unreachable!(),
            };

            println!("{} > {:?}", path_seg, nodes);

            let encoded_path = encode_path(nodes, &rules, &mut statistics, &mut statistics2, path_id, offset);
            result.insert(path_seg.to_string(), encoded_path);

            path_id += 1;

        }
        buf.clear();
    }
    for (key, usage) in statistics {
        if key[0] < offset && key[1] < offset {
        if !rules[&key].colors.equal_set(&usage) || rules[&key].colors.len() != statistics2[&key] {
            log::warn!("Used {} ({}) - {} ({}) exactly {} times vs supposed {}, difference: {} | {}", key[0], key[0] >= offset, key[1], key[1] >= offset, statistics2[&key], rules[&key].colors.len(), usage.difference(&rules[&key].colors), rules[&key].colors.difference(&usage));
            if rules[&key].colors.len() < statistics2[&key] {
                log::warn!("usage: {}", usage);
            }
        }
        }
    }
    result
}

pub fn encode_path(nodes: Vec<NodeId>, rules: &HashMap<Vec<NodeId>, Rule>, statistics: &mut HashMap<Vec<NodeId>, ColorSet>, statistics2: &mut HashMap<Vec<NodeId>, usize>, path_id: u64, offset: NodeId) -> Vec<NodeId> {
    let mut stack: Vec<NodeId> = Vec::new();
    let mut multiplicity_stack: Vec<(usize, usize)> = Vec::new();

    let mut nodes_visited_prev: HashSet<NodeId> = HashSet::new();
    let mut prev_counter = 0;
    let mut nodes_visited_curr: HashSet<NodeId> = HashSet::new();
    let mut curr_counter = 0;

    log::info!("Encoding path {}", path_id);

    for (idx, current_node) in nodes.into_iter().enumerate() {
        match stack.pop() {
            None => {
                stack.push(current_node);
                nodes_visited_prev.insert(current_node.get_forward());
            },
            Some(prev_node) => {
                if nodes_visited_curr.contains(&current_node.get_forward()) {
                    curr_counter += 1;
                    nodes_visited_curr.clear();
                }
                nodes_visited_curr.insert(current_node.get_forward());

                let mut has_used_rule = false;

                let digram = canonize(prev_node, current_node);
                let (first_counter, second_counter) = (prev_counter, curr_counter);
                if let Some(rule) = rules.get(&vec![digram.0, digram.1]) {
                    if digram.0 == NodeId::new(5, 0) || digram.1 == NodeId::new(5, 1) {
                        log::error!("Using rule {} {}:{} twice at {} | {} - {}", rule, prev_node, current_node, idx, first_counter, second_counter);
                    }
                    if rule.colors.contains(path_id, first_counter, second_counter) {
                        if statistics[&vec![digram.0, digram.1]].contains(path_id, first_counter, second_counter) {
                            log::error!("Using rule {} twice at {}", rule, idx);
                        }
                        if digram.0 == prev_node {
                            stack.push(rule.left);
                        } else {
                            stack.push(rule.left.flip());
                        }
                        statistics.get_mut(&vec![digram.0, digram.1]).expect("statistics contains rule").insert(path_id, first_counter, second_counter);
                        *statistics2.get_mut(&vec![digram.0, digram.1]).expect("statistics contains rule") += 1;
                        has_used_rule = true;
                    } else {
                        // log::debug!("rule: {}, {},{} | color: {},{},{} | colors: {}", rule, prev_node, current_node, path_id, first_counter, second_counter, rule.colors);
                    }
                }

                if !has_used_rule {
                    stack.push(prev_node);
                    stack.push(current_node);
                    multiplicity_stack.push((first_counter, second_counter));
                } else if stack.len() >= 2 {
                    let mut change = true;
                    while change && stack.len() >= 2 {
                        let second = stack.pop().expect("stack has at least one node");
                        let first = stack.pop().expect("stack has at least two nodes");
                        let mult = multiplicity_stack.pop().expect("multiplicity stack has at least one entry");
                        let digram = canonize(first, second);
                        change = false;
                        if let Some(rule) = rules.get(&vec![digram.0, digram.1]) {
                            if rule.colors.contains(path_id, mult.0, mult.1) {
                                if digram.0 == first {
                                    stack.push(rule.left);
                                } else {
                                    stack.push(rule.left.flip());
                                }
                                statistics.get_mut(&vec![digram.0, digram.1]).expect("statistics contains rule").insert(path_id, mult.0, mult.1);
                                *statistics2.get_mut(&vec![digram.0, digram.1]).expect("statistics contains rule") += 1;
                                change = true;
                            }
                        }
                        if !change {
                            stack.push(first);
                            stack.push(second);
                            multiplicity_stack.push(mult);
                        }
                    }
                }

                // Do this last, to set prev only at the end and only once
                if nodes_visited_prev.contains(&current_node.get_forward()) {
                    prev_counter += 1;
                    nodes_visited_prev.clear()
                }
                nodes_visited_prev.insert(current_node.get_forward());
            },
        }
    }
    stack
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
            parents.entry(u).and_modify(|e| *e = None).or_insert(Some(non_terminal));
        }
        if v >= offset {
            parents.entry(u).and_modify(|e| *e = None).or_insert(Some(non_terminal));
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

        if u == NodeId::new(26, 0) || u == NodeId::new(26, 1){
            log::error!("We are at the critical point - {} | {}", u, v);
        }
        // if v == NodeId::new(8, 0) {
        //     log::error!("We are at the critical point - {} | {}", u, v);
        // }

        for n in &neighbors[v.get_idx()] {
            let n = *n;
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
            // if v == NodeId::new(5, 0) {
            //     log::error!("vn - {} | {} | {} - {} - {}", v, n, digrams.get_priority(&canonize(v, n)).unwrap(), qn_set, vn_set);
            //     log::error!("qn - {} - {} - {}", digrams.get_priority(&canonize(v, n)).unwrap(), uv_color_set, qn_set);
            // }
            // if v == NodeId::new(8, 0) {
            //     log::error!("vn - {} | {} | {} - {} - {}", v, n, digrams.get_priority(&canonize(v, n)).unwrap(), qn_set, vn_set);
            //     log::error!("qn - {} - {} - {}", digrams.get_priority(&canonize(v, n)).unwrap(), uv_color_set, qn_set);
            // }
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
            if u == NodeId::new(26, 1) || u == NodeId::new(26, 0) {
                log::error!("nu - {} | {} | {} - {}", n, u, digrams.get_priority(&canonize(n, u)).unwrap(), nq_set);
                log::error!("nq - {} - {} - {}", digrams.get_priority(&canonize(n, u)).unwrap(), uv_color_set, nq_set);
            }
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
        }
        neighbors_to_insert.into_iter().for_each(|(key, value)| {
            neighbors[key.get_idx()].insert(value);
        });
        neighbors_to_remove.into_iter().for_each(|(key, value)| {
            neighbors[key.get_idx()].remove(&value);
        });

        rules.insert(non_terminal, Rule { left: non_terminal, right: vec![u, v], colors: uv_color_set });
    }
    log::info!("Built {} rules, more than 0 parents: {}, exactly 1 parent: {}", rules.len(), parents.len(), parents.iter().map(|(_, x)| if x.is_some() { 1 } else { 0 }).sum::<u64>());

    let rules: Rules = rules.into_iter().filter(|(_k, v)| !v.colors.is_empty() || v.colors.len() == 1).collect();
    log::info!("After removal: {}", rules.len());

    let rules = merge_rules(rules, parents, offset);

    (rules, offset)
}

fn reverse_rule(right: &Vec<NodeId>) -> Vec<NodeId> {
    right.iter().copied().rev().map(|x| x.flip()).collect()
}

fn merge_rules(mut rules: Rules, parents: HashMap<NodeId, Option<NodeId>>, offset: NodeId) -> Rules {
    log::info!("Merging rules");
    let mut mergeables: VecDeque<NodeId> = VecDeque::new();
    let mut counter = 0;
    let mut counter_diff_colors = 0;
    let mut counter_not_one_parent = 0;
    for (key, rule) in rules.iter() {
        if rule.right.iter().all(|&x| x < offset) {
            mergeables.push_back(*key);
        }
    }
    log::debug!("Mergeables: {}", mergeables.len());
    while !mergeables.is_empty() {
        let rule_to_collapse = mergeables.pop_front().expect("mergeables should contain at least one element");
        if let Some(parent) = parents.get(&rule_to_collapse).cloned().flatten() {
            if rules[&rule_to_collapse].colors.equal_set(&rules[&parent].colors) {
                let mut idx = 0;
                let mut reverse = false;
                for (index, node) in rules[&parent].right.iter().enumerate() {
                    if node.get_forward() == rule_to_collapse {
                        idx = index;
                        reverse = node.get_orientation() == 1;
                    }
                }
                let rule_to_insert = if reverse {
                    reverse_rule(&rules[&rule_to_collapse].right)
                } else {
                    rules[&rule_to_collapse].right.clone()
                };
                rules.get_mut(&parent).expect("rules has parent").right.splice(idx..idx + 1, rule_to_insert);
                rules.remove(&rule_to_collapse).expect("rules contains rule to remove");
                if rules[&parent].right.iter().all(|&x| x < offset) {
                    mergeables.push_back(parent);
                }
                counter += 1;
            } else {
                counter_diff_colors += 1;
            }
        } else {
            counter_not_one_parent += 1;
        }
    }
    log::info!("Merged {} rules, {} times the color set was different, {} times there was not one parent", counter, counter_diff_colors, counter_not_one_parent);
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
        println!("{} - {}", rule, rule.colors);
    }
    for result in encoded_paths {
        println!("{:?}", result);
    }
}
