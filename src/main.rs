mod color_set;
mod node_id;
mod parser;
mod path_segment;

use clap::Parser;
use color_set::ColorSet;
use node_id::{NodeId, RawNodeId};
use parser::{
    bufreader_from_compressed_gfa, canonize, get_nodes_path, get_nodes_walk, parse_gfa_paths_walks,
    parse_node_ids, parse_path_identifier, parse_walk_identifier,
};
use priority_queue::PriorityQueue;
use std::collections::{HashMap, HashSet};
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
        write!(f, "{} -> {:?}", self.left, self.right,)
    }
}

impl fmt::Display for Rule {
    // This trait requires `fmt` with this exact signature.
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{} -> {:?}", self.left, self.right,)
    }
}

pub fn encode_paths(
    file: &str,
    rules: &Rules,
    offset: NodeId,
    node_ids_by_name: &HashMap<Vec<u8>, RawNodeId>,
) -> HashMap<String, Vec<NodeId>> {
    let mut result: HashMap<String, Vec<NodeId>> = HashMap::new();
    let mut rules_by_digram: HashMap<(NodeId, NodeId), Vec<Rule>> = HashMap::new();
    for (_idx, rule) in rules.iter() {
        rules_by_digram
            .entry(canonize(rule.right[0], rule.right[1]))
            .and_modify(|e| e.push(rule.clone()))
            .or_insert(vec![rule.clone()]);
    }
    let mut rules_by_end_digram: HashMap<(NodeId, NodeId), Vec<&Rule>> = HashMap::new();
    for (_idx, rule) in rules.iter() {
        rules_by_end_digram
            .entry(canonize(rule.right[0], rule.right[1]))
            .and_modify(|e| e.push(rule))
            .or_insert(vec![rule]);
    }
    let mut data = bufreader_from_compressed_gfa(file);
    let mut statistics: HashMap<NodeId, usize> =
        rules.iter().map(|(k, _v)| (k.clone(), 0)).collect();
    let mut path_id = 0;

    log::error!("RBD: {:?}", rules_by_digram);

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

            let encoded_path = encode_path(
                nodes,
                &rules_by_digram,
                &rules_by_end_digram,
                &mut statistics,
                path_id,
                offset,
            );
            result.insert(path_seg.to_string(), encoded_path);

            path_id += 1;
        }
        buf.clear();
    }
    for (key, usage) in statistics {
        if rules[&key].colors.len() != usage {
            log::warn!(
                "Used {:?} exactly {} times vs supposed {}",
                key,
                usage,
                rules[&key].colors.len()
            );
            // break;
        }
    }
    result
}

pub fn encode_path(
    nodes: Vec<NodeId>,
    rules: &HashMap<(NodeId, NodeId), Vec<Rule>>,
    rules_by_end: &HashMap<(NodeId, NodeId), Vec<&Rule>>,
    statistics: &mut HashMap<NodeId, usize>,
    path_id: u64,
    _offset: NodeId,
) -> Vec<NodeId> {
    let mut stack: Vec<NodeId> = Vec::new();
    let mut multiplicity_stack: Vec<(usize, usize)> = Vec::new();

    let mut nodes_visited_prev: HashSet<NodeId> = HashSet::new();
    let mut prev_counter = 0;
    let mut nodes_visited_curr: HashSet<NodeId> = HashSet::new();
    let mut curr_counter = 0;

    log::info!("Encoding path {}", path_id);

    let mut idx = 0;
    while idx < nodes.len() {
        let current_node = nodes[idx];
        log::error!("Working on node {}, idx {}", current_node, idx);
        match stack.pop() {
            None => {
                stack.push(current_node);
                nodes_visited_prev.insert(current_node.get_forward());
            }
            Some(prev_node) => {
                if nodes_visited_curr.contains(&current_node.get_forward()) {
                    log::error!("= Resetting outside curr =");
                    curr_counter += 1;
                    nodes_visited_curr.clear();
                }
                nodes_visited_curr.insert(current_node.get_forward());

                let (first_counter, second_counter) = (prev_counter, curr_counter);
                let has_used_rule = apply_rule(
                    &mut idx,
                    (prev_node, current_node),
                    (first_counter, second_counter),
                    &mut stack,
                    &mut multiplicity_stack,
                    &nodes,
                    statistics,
                    rules,
                    rules_by_end,
                    path_id,
                    &mut nodes_visited_prev,
                    &mut nodes_visited_curr,
                    &mut prev_counter,
                    &mut curr_counter,
                );

                if !has_used_rule {
                    log::error!("Has not used rule for {} {}", prev_node, current_node);
                    stack.push(prev_node);
                    stack.push(current_node);
                    multiplicity_stack.push((first_counter, second_counter));
                } else if stack.len() >= 2 {
                    let mut change = true;
                    while change && stack.len() >= 2 {
                        let second = stack.pop().expect("stack has at least one node");
                        let first = stack.pop().expect("stack has at least two nodes");
                        let mult = multiplicity_stack
                            .pop()
                            .expect("multiplicity stack has at least one entry");
                        change = apply_rule(
                            &mut idx,
                            (first, second),
                            mult,
                            &mut stack,
                            &mut multiplicity_stack,
                            &nodes,
                            statistics,
                            rules,
                            rules_by_end,
                            path_id,
                            &mut nodes_visited_prev,
                            &mut nodes_visited_curr,
                            &mut prev_counter,
                            &mut curr_counter,
                        );
                        if !change {
                            stack.push(first);
                            stack.push(second);
                            multiplicity_stack.push(mult);
                        }
                    }
                }

                // Do this last, to set prev only at the end and only once
                if nodes_visited_prev.contains(&current_node.get_forward()) {
                    log::error!("= Resetting outside prev =");
                    prev_counter += 1;
                    nodes_visited_prev.clear()
                }
                nodes_visited_prev.insert(current_node.get_forward());
            }
        }
        log::error!("stack: {:?}, idx: {}", stack, idx);
        idx += 1;
    }
    stack
}

fn apply_rule(
    idx: &mut usize,
    digram: (NodeId, NodeId),
    multiplicity: (usize, usize),
    stack: &mut Vec<NodeId>,
    multiplicity_stack: &mut Vec<(usize, usize)>,
    nodes: &Vec<NodeId>,
    statistics: &mut HashMap<NodeId, usize>,
    rules: &HashMap<(NodeId, NodeId), Vec<Rule>>,
    rules_by_end: &HashMap<(NodeId, NodeId), Vec<&Rule>>,
    path_id: u64,
    nodes_visited_prev: &mut HashSet<NodeId>,
    nodes_visited_curr: &mut HashSet<NodeId>,
    prev_counter: &mut usize,
    curr_counter: &mut usize,
) -> bool {
    let canonized_digram = canonize(digram.0, digram.1);
    log::error!("Canonized: {:?}", canonized_digram);
    let mut has_used_rule = false;
    for (is_end, rule) in rules
        .get(&canonized_digram)
        .into_iter()
        .flatten()
        .map(|r| (false, r))
        .chain(
            rules_by_end
                .get(&canonized_digram)
                .into_iter()
                .flatten()
                .map(|r| (true, *r)),
        )
    {
        log::error!("\t-- Trying to work on rule: {}, with colors {}, mult: {:?} --", rule, rule.colors, multiplicity);
        if rule
            .colors
            .contains(path_id, multiplicity.0, multiplicity.1)
        {
            log::error!("\t== Working on rule: {} ==", rule);
            if (!is_end && digram.0 == rule.right[0]) || (is_end && digram.0 == rule.right[rule.right.len() - 1].flip()) {
                // Exact match, continue forward
                log::error!(
                    "{:?} - forward match | nodes_len {} | idx {} | rule_len {}",
                    digram,
                    nodes.len(),
                    idx,
                    rule.right.len()
                );
                if rule.right.len() > 2 && nodes.len() > *idx + rule.right.len() - 2 {
                    log::error!(
                        "\tRemaining rule {:?} - {:?}",
                        &rule.right[2..],
                        &nodes[*idx + 1..*idx + 1 + rule.right.len() - 2]
                    );
                }
                if (rule.right.len() > 2
                    && !is_end
                    && nodes.len() - *idx > 1
                    && rule.right[2..] == nodes[*idx + 1..*idx + 1 + rule.right.len() - 2])
                    || (rule.right.len() > 2
                        && is_end
                        && nodes.len() - *idx > 1
                        && rule.right[..rule.right.len() - 2] == nodes[*idx + 1..*idx + 1 + rule.right.len() - 2])
                    || (rule.right.len() == 2)
                {
                    // Forward match
                    if !is_end {
                        stack.push(rule.left.clone());
                    } else {
                        stack.push(rule.left.flip());
                    }
                    log::error!("\tInner forward: stack: {:?}", stack);
                    *statistics
                        .get_mut(&rule.left)
                        .expect("statistics contains rule") += 1;
                    has_used_rule = true;

                    let new_idx = *idx + rule.right.len() - 2;
                    for i in *idx + 1..=new_idx {
                        if nodes_visited_curr.contains(&nodes[i].get_forward()) {
                            log::error!("= Resetting inside curr =");
                            *curr_counter += 1;
                            nodes_visited_curr.clear();
                        }
                        nodes_visited_curr.insert(nodes[i].get_forward());
                        if nodes_visited_prev.contains(&nodes[i].get_forward()) {
                            log::error!("= Resetting inside prev =");
                            *prev_counter += 1;
                            nodes_visited_prev.clear()
                        }
                        nodes_visited_prev.insert(nodes[i].get_forward());
                    }
                    *idx = new_idx;
                    break;
                }
            } else if (!is_end && digram.0 == rule.right[1].flip()) || (is_end && digram.0 == rule.right[rule.right.len() - 2]) {
                // Reverse match continue backward
                let rule_remaining_len = rule.right.len() - 2;
                log::error!(
                    "{:?} - backward match, is_end: {},  rule_remaining_len {}, stack {} -> {} {}",
                    digram,
                    is_end,
                    rule_remaining_len,
                    stack.len(),
                    rule_remaining_len > 0,
                    stack.len() >= 1 + rule_remaining_len
                );
                if rule_remaining_len > 0 && stack.len() >= rule_remaining_len {
                    log::error!(
                        "Remaining rule {:?} - {:?}",
                        &rule.right[2..],
                        reverse_rule(&stack[stack.len() - rule_remaining_len..].to_vec())
                    );
                }
                if (rule_remaining_len == 0)
                    || (rule_remaining_len > 0
                        && !is_end
                        && stack.len() >= rule_remaining_len
                        && rule.right[2..]
                            == reverse_rule(&stack[stack.len() - rule_remaining_len..].to_vec()))
                    || (rule_remaining_len > 0
                        && is_end
                        && stack.len() >= rule_remaining_len
                        && rule.right[..rule_remaining_len]
                            == stack[stack.len() - rule_remaining_len..])
                {
                    log::error!("\tInner backward");
                    for _ in 0..rule_remaining_len {
                        stack.pop();
                        multiplicity_stack.pop();
                    }

                    if !is_end {
                        stack.push(rule.left.flip());
                    } else {
                        stack.push(rule.left);
                    }
                    *statistics
                        .get_mut(&rule.left)
                        .expect("statistics contains rule") += 1;
                    has_used_rule = true;
                    break;
                }
            }
        } else {
            log::error!(
                "\tNo color match: color: {}, {} | {} | {}",
                rule.colors,
                path_id,
                multiplicity.0,
                multiplicity.1
            );
        }
    }
    has_used_rule
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

        for n in &neighbors[v.get_idx()] {
            let n = *n;
            if n == v && u == v {
                // log::debug!("Self loop case: {} -> {} - {} - {}", non_terminal, u, v, n);
                neighbors_to_remove.push((v, n));
                neighbors_to_remove.push((n.flip(), v.flip()));
                continue;
            }
            if digrams.get_priority(&canonize(v, n)).is_none() {
                log::error!("Missing digram: {}-{}, for rule: {}-{}", v, n, u, v);
            }
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
            let mut empty_vn_set = false;
            if vn_set.is_empty() {
                neighbors_to_remove.push((v, n));
                neighbors_to_remove.push((n.flip(), v.flip()));
                empty_vn_set = true;
            }

            if n == u && u != v {
                // log::debug!("Loop over rule case: {} -> {} - {} - {}", non_terminal, u, v, n);
                digrams.push(canonize(non_terminal, non_terminal), qn_set);
                digrams.change_priority(&canonize(v, n), vn_set);

                // <n += <(u>v>)
                neighbors_to_insert.push((non_terminal.flip(), non_terminal.flip()));
                neighbors_to_insert.push((non_terminal, non_terminal));
            } else {
                digrams.push(canonize(non_terminal, n), qn_set);
                digrams.change_priority(&canonize(v, n), vn_set);

                // <n += <(u>v>)
                neighbors_to_insert.push((n.flip(), non_terminal.flip()));
                neighbors_to_insert.push((non_terminal, n));
            }

            if empty_vn_set {
                digrams.remove(&canonize(v, n));
            }
        }
        for n in neighbors.get(u.flip().get_idx()).unwrap() {
            let n = n.flip();
            if n == v {
                continue;
            }
            let nq_set = digrams
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
        }
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
        println!("Q\t{} | {}", rule, rule.colors);
    }
    for (path_name, path) in encoded_paths {
        println!("Z\t{}\t{:?}", path_name, path);
    }
}
