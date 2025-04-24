mod color_set;
mod compressor;
mod node_id;
mod parser;
mod path_segment;

use clap::{Parser, Subcommand};
use color_set::OrdColorSet;
use compressor::encode_paths2;
use node_id::{NodeId, RawNodeId};
use parser::{
    bufreader_from_compressed_gfa, canonize, parse_compressed_lines, parse_gfa_paths_walks,
    parse_node_ids,
};
use path_segment::PathSegment;
use priority_queue::PriorityQueue;
use std::collections::{HashMap, HashSet};
use std::fmt;
use std::io::BufRead;
use std::mem;

const MAX_OCCURENCES: usize = 2;

type NeighborList = Vec<HashSet<NodeId>>;
type Digrams = PriorityQueue<(NodeId, NodeId), OrdColorSet>;
type Rules = HashMap<NodeId, Rule>;

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

pub fn build_qlines(
    neighbors: &mut NeighborList,
    digrams: &mut Digrams,
) -> (Rules, NodeId, HashMap<NodeId, Vec<NodeId>>) {
    log::info!("Building qlines for {} digrams", digrams.len());
    let offset = NodeId::from_raw(neighbors.len() as u64);
    let mut rules: Rules = HashMap::new();
    let mut parents: HashMap<NodeId, Vec<NodeId>> = HashMap::new();

    let mut current_max_node_id = offset;

    while digrams.peek().expect("At least one digram").1.len() >= MAX_OCCURENCES {
        let ((u, v), uv_color_set) = digrams.pop().expect("At least one digram");
        let non_terminal: NodeId = current_max_node_id;
        current_max_node_id += 2;

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
        let should_print = false;

        for n in neighbors.get(u.flip().get_idx()).unwrap() {
            let n = n.flip();
            // println!("nu's n: {}", n);

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
                    println!(
                        "non_terminal: {}, nq_set: {:?}, mut: {:?}",
                        non_terminal, nq_set, mutation_outgoing
                    );
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
            let self_sets = uv_color_set.colors.sectionize(&mut mutation_outgoing);
            // log::error!("Self loop contents: qq: {:?}, qv: {:?}, uv: {:?}", self_sets.0, self_sets.1, self_sets.2);
            Some(self_sets)
        } else {
            None
        };

        for n in &neighbors[v.get_idx()] {
            let n = *n;
            // println!("vn's n: {}", n);

            if n == u && u != v {
                // TODO: handle this case
            } else if n == v && u == v {
                let (mut qq_set, mut qv_set, uv_temp) =
                    mem::take(&mut self_sets).expect("self sets should have been set");
                // println!("qq_set: {:?}", qq_set);
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
                let (new_qn_set, vn_set_addition) =
                    qn_set.self_vy_intersection(&mutation_outgoing, is_vn_flipped, is_qn_flipped);
                qn_set = new_qn_set;
                new_vn_set.add_addition(vn_set_addition);
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

        if u >= offset {
            parents
                .entry(u.get_forward())
                .and_modify(|e| e.push(non_terminal))
                .or_insert(vec![non_terminal]);
        }
        if v >= offset && u != v {
            parents
                .entry(v.get_forward())
                .and_modify(|e| e.push(non_terminal))
                .or_insert(vec![non_terminal]);
        }

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
            // println!("Self sets: {:?}", self_sets);
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

    let rules: Rules = rules
        .into_iter()
        // .filter(|(_k, v)| !v.colors.colors.is_empty() || v.colors.len() == 1)
        .collect();
    log::info!("Built {} rules", rules.len(),);
    (rules, offset, parents)
}

fn is_edge_flipped(a: NodeId, b: NodeId) -> bool {
    canonize(a, b).0 == b.flip()
}

fn insert_edge(neighbors_to_insert: &mut Vec<(NodeId, NodeId)>, a: NodeId, b: NodeId) {
    neighbors_to_insert.push((a, b));
    neighbors_to_insert.push((b.flip(), a.flip()));
}

fn reverse_rule(right: &Vec<NodeId>) -> Vec<NodeId> {
    right.iter().copied().rev().map(|x| x.flip()).collect()
}

fn is_mergeable(
    rule: NodeId,
    rules: &Rules,
    parents: &HashMap<NodeId, Vec<NodeId>>,
    threshold: usize,
) -> bool {
    let rule_uses = parents.get(&rule).map(|p| p.len()).unwrap_or_default();
    let all_uses = rules[&rule].colors.colors.len();
    let path_uses = all_uses
        - parents
            .get(&rule)
            .map(|rule_parents| {
                rule_parents
                    .iter()
                    .map(|parent| {
                        if !rules.contains_key(parent) {
                            log::error!(
                                "Missing rules for parent: {} (parent of {}, content: {:?})",
                                parent,
                                rule,
                                rules[&rule].right
                            );
                        }
                        rules[parent].colors.colors.len()
                    })
                    .sum::<usize>()
            })
            .unwrap_or_default();
    rule_uses + path_uses < threshold
}

fn replace_rule(text: &mut Vec<NodeId>, rule: NodeId, replacement: Vec<NodeId>) {
    let mut idx = 0;
    let mut reverse = false;
    let mut no_hit = true;
    for (index, node) in text.iter().enumerate() {
        if node.get_forward() == rule {
            idx = index;
            reverse = node.get_orientation() == 1;
            no_hit = false;
        }
    }
    if no_hit {
        log::warn!("Cannot replace rule: {}", rule);
        return;
    }
    let rule_to_insert = if reverse {
        reverse_rule(&replacement)
    } else {
        replacement
    };
    text.splice(idx..idx + 1, rule_to_insert);
}

fn get_correct_path_ids(
    parents: &HashMap<NodeId, Vec<NodeId>>,
    rules: &Rules,
    rule: NodeId,
) -> Vec<u64> {
    let mut all_path_ids = rules[&rule].colors.colors.get_path_counts();
    for parent in parents.get(&rule).unwrap_or(&Vec::new()).iter() {
        let parent_path_ids = rules[parent].colors.colors.get_path_counts();
        for (path_id, count) in parent_path_ids {
            if all_path_ids[&path_id] >= count {
                *all_path_ids
                    .get_mut(&path_id)
                    .expect("Child contains all paths of parent") -= count;
            } else {
                log::error!("Getting a negative amount of paths");
            }
        }
    }
    all_path_ids
        .into_iter()
        .filter_map(|(path_id, count)| if count > 0 { Some(path_id) } else { None })
        .collect::<Vec<u64>>()
}

fn simplify_rules(
    mut rules: Rules,
    parents: &mut HashMap<NodeId, Vec<NodeId>>,
    encoded_paths: &mut HashMap<PathSegment, Vec<NodeId>>,
    path_id_to_path_segment: &HashMap<u64, PathSegment>,
    threshold: usize,
) -> Rules {
    log::info!("Simplifying rules");

    let mut potentially_mergeables: Vec<NodeId> = rules.keys().copied().collect();
    let mut counter = 0;

    while !potentially_mergeables.is_empty() {
        let rule = potentially_mergeables
            .pop()
            .expect("potentially_mergeables has at least one item");
        if !is_mergeable(rule, &rules, &parents, threshold) {
            continue;
        }
        log::debug!("Simplifying rule: {}", rule);
        for parent in parents.get(&rule).unwrap_or(&Vec::new()).iter() {
            let replacement = rules[&rule].right.clone();
            replace_rule(
                &mut rules.get_mut(parent).expect("Rules contain parent").right,
                rule,
                replacement,
            );
        }
        for path_id in get_correct_path_ids(&parents, &rules, rule) {
            let path_name = &path_id_to_path_segment[&path_id];
            let replacement = rules[&rule].right.clone();
            replace_rule(
                encoded_paths
                    .get_mut(&path_name)
                    .expect("Encoded paths contains path"),
                rule,
                replacement,
            );
        }
        for child in rules[&rule].right.iter() {
            let child = child.get_forward();
            if rules.contains_key(&child) {
                parents
                    .get_mut(&child)
                    .expect("Parents of child contain node")
                    .retain(|el| *el != rule);
                let rule_parents = parents.get(&rule).cloned().unwrap_or_default();
                parents
                    .get_mut(&child)
                    .expect("Parents of child contain node")
                    .extend(rule_parents);
                // if !potentially_mergeables.contains(&child) {
                //     potentially_mergeables.push(child);
                // }
            }
        }
        rules.remove(&rule);
        counter += 1;
    }
    log::info!("Removed {} rules during simplification", counter);
    rules
}

fn check_rule_usability(
    offset: NodeId,
    encoded_paths: &HashMap<PathSegment, Vec<NodeId>>,
    rules: &Rules,
    parents: &HashMap<NodeId, Vec<NodeId>>,
) -> () {
    let mut rule_usage_path: HashMap<NodeId, usize> = HashMap::new();
    let mut rule_usage_rule: HashMap<NodeId, usize> = HashMap::new();
    let mut everything_ok = true;
    for (_path_name, path) in encoded_paths {
        for node in path {
            if *node >= offset {
                rule_usage_path
                    .entry(node.get_forward())
                    .and_modify(|e| *e += 1)
                    .or_insert(1);
                if !rules.contains_key(&node.get_forward()) {
                    log::error!(
                        "Rule {} was mistakenly deleted (used in path)",
                        node.get_forward()
                    );
                    everything_ok = false;
                }
            }
        }
    }
    for (_rule_name, rule) in rules {
        for node in &rule.right {
            if *node >= offset {
                rule_usage_rule
                    .entry(node.get_forward())
                    .and_modify(|e| *e += 1)
                    .or_insert(1);
                if !rules.contains_key(&node.get_forward()) {
                    log::error!(
                        "Rule {} was mistakenly deleted (used in other rule)",
                        node.get_forward()
                    );
                    everything_ok = false;
                }
            }
        }
    }
    let used_rules: HashSet<NodeId> = rule_usage_path
        .keys()
        .chain(rule_usage_rule.keys())
        .copied()
        .collect();
    for rule in &used_rules {
        if rule_usage_path.get(rule).copied().unwrap_or_default()
            + rule_usage_rule.get(rule).copied().unwrap_or_default()
            == 1
        {
            let path_appearance = rule_usage_path.get(rule).copied().unwrap_or_default() == 1;
            if parents.contains_key(rule) {
                log::warn!(
                    "Rule {} can be removed, but wasn't (path: {}, parents: {:?}, colors: {:?})",
                    rule,
                    path_appearance,
                    parents[rule],
                    rules[rule].colors
                );
            } else {
                log::warn!(
                    "Rule {} can be removed, but wasn't (path: {}, parents: zero, colors: {:?})",
                    rule,
                    path_appearance,
                    rules[rule].colors
                );
            }
            everything_ok = false;
        }
    }

    if everything_ok {
        log::info!("Rule usability is optimal");
    }
}

fn write_compressed_lines(gfa_file: &str, rules: &Rules, encoded_paths: &HashMap<PathSegment, Vec<NodeId>>, node_ids_by_name: HashMap<Vec<u8>, RawNodeId>) {
    log::info!("Writing output");
    let mut buf = vec![];
    let mut data = bufreader_from_compressed_gfa(gfa_file);

    let node_names_by_id: HashMap<RawNodeId, Vec<u8>> = node_ids_by_name.into_iter().map(|(k, v)| (v, k)).collect();

    // Print all other lines
    while data.read_until(b'\n', &mut buf).unwrap_or(0) > 0 {
        if buf[0] == b'P' || buf[0] == b'W' {
            buf.clear();
            continue;
        } else {
            unsafe {
                print!("{}", std::str::from_utf8_unchecked(&buf));
            }
            buf.clear();
        }
    }

    // Print all Q-lines
    for rule in rules.values() {
        unsafe {
            print!("Q\t{}\t", std::str::from_utf8_unchecked(&node_names_by_id[&rule.left.get_id()]));
        }
        for node in &rule.right {
            unsafe {
                print!("{}{}", if node.is_forward() { ">" } else { "<" }, std::str::from_utf8_unchecked(&node_names_by_id[&node.get_id()]));
            }
        }
        println!();
    }

    // Print all Z-lines
    for (path_name, path) in encoded_paths {
        print!("Z\t{}\t", path_name.to_walk_string());
        for node in path {
            unsafe {
                print!("{}{}", if node.is_forward() { ">" } else { "<" }, std::str::from_utf8_unchecked(&node_names_by_id[&node.get_id()]));
            }
        }
        println!();
    }
}

fn write_decompressed_lines(gfa_file: &str, decompressed_paths: HashMap<PathSegment, Vec<NodeId>>, node_ids_by_name: HashMap<Vec<u8>, RawNodeId>, use_p_lines: bool) {
    log::info!("Writing output");
    let mut buf = vec![];
    let mut data = bufreader_from_compressed_gfa(gfa_file);

    let node_names_by_id: HashMap<RawNodeId, Vec<u8>> = node_ids_by_name.into_iter().map(|(k, v)| (v, k)).collect();

    // Print all other lines
    while data.read_until(b'\n', &mut buf).unwrap_or(0) > 0 {
        if buf[0] == b'Q' || buf[0] == b'Z' {
            buf.clear();
            continue;
        } else {
            unsafe {
                print!("{}", std::str::from_utf8_unchecked(&buf));
            }
            buf.clear();
        }
    }

    // Print all Z-lines
    for (path_name, path) in decompressed_paths {
        if !use_p_lines {
            print!("W\t{}\t", path_name.to_walk_string());
            for node in path {
                unsafe {
                    print!("{}{}", if node.is_forward() { ">" } else { "<" }, std::str::from_utf8_unchecked(&node_names_by_id[&node.get_id()]));
                }
            }
            println!();
        } else {
            print!("P\t{}\t", path_name);
            let mut iter = path.iter();
            if let Some(node) = iter.next() {
                unsafe {
                    print!("{}{}", std::str::from_utf8_unchecked(&node_names_by_id[&node.get_id()]), if node.is_forward() { "+" } else { "-" });
                }
            }
            for node in iter {
                unsafe {
                    print!(",{}{}", std::str::from_utf8_unchecked(&node_names_by_id[&node.get_id()]), if node.is_forward() { "+" } else { "-" });
                }
            }
            println!();
        }
    }
}

fn get_sequence(node: NodeId, rules: &HashMap<NodeId, Vec<NodeId>>) -> Vec<NodeId> {
    if let Some(rule) = rules.get(&node.get_forward()) {
        if node.is_forward() {
            rule.iter().flat_map(|n| get_sequence(*n, rules)).collect()
        } else {
            reverse_rule(&rule.iter().flat_map(|n| get_sequence(*n, rules)).collect())
        }
    } else {
        vec![node]
    }
}

fn decompress(rules: HashMap<NodeId, Vec<NodeId>>, paths: HashMap<PathSegment, Vec<NodeId>>) -> HashMap<PathSegment, Vec<NodeId>> {
    let mut decompressed_paths: HashMap<PathSegment, Vec<NodeId>> = HashMap::new();
    for (path_name, path) in paths {
        let mut new_sequence = Vec::new();
        for node in path {
            new_sequence.extend(get_sequence(node, &rules));
        }
        decompressed_paths.insert(path_name, new_sequence);
    }
    decompressed_paths
}

fn get_non_terminal_node_names(rules: &HashMap<NodeId, Rule>, prefix: String) -> HashMap<Vec<u8>, RawNodeId> {
    let mut result = HashMap::new();
    let mut counter = 1;
    for non_terminal in rules.keys() {
        let name = format!("{}{}", prefix, counter);
        result.insert(name.into_bytes(), non_terminal.get_id());
        counter += 1;
    }
    result
}

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Debug, Subcommand)]
enum Commands {
    /// Compresses a GFA file
    #[command(arg_required_else_help = true)]
    Compress {
        /// Input GFA file
        #[arg(required = true)]
        file: String,

        /// Minimum number a sequence has to appear to have its own rule
        #[arg(short, default_value_t = 2)]
        k: usize,

        /// Prefix that is used for non-terminal identifiers (i.e. Q-lines)
        #[arg(short, long, default_value = "Q")]
        prefix: String,
    },

    /// Decompresses a GFA file (by default paths will be decompressed as W-lines)
    #[command(arg_required_else_help = true)]
    Decompress {
        /// Input GFA file
        #[arg(required = true)]
        file: String,

        #[arg(short = 'p', long)]
        use_p_lines: bool,
    },
}

fn main() {
    env_logger::init();
    let args = Cli::parse();

    match args.command {
        Commands::Compress { file, k, prefix } => {
            log::info!("Compressing file {} with k = {}", file, k);
            let mut node_ids_by_name = parse_node_ids(&file, false);
            let (mut neighbors, mut digrams, path_id_to_path_segment) =
                parse_gfa_paths_walks(&file, &node_ids_by_name);
            let (rules, offset, parents) = build_qlines(&mut neighbors, &mut digrams);
            let mut encoded_paths = encode_paths2(&file, &rules, offset, &node_ids_by_name);
            let mut prules = rules.iter().collect::<Vec<_>>();
            prules.sort_by_key(|r| r.0.get_idx());
            let mut safe_parents = parents.clone();
            let rules = simplify_rules(
                rules,
                &mut safe_parents,
                &mut encoded_paths,
                &path_id_to_path_segment,
                k,
            );
            check_rule_usability(offset, &encoded_paths, &rules, &parents);
            let non_terminal_node_names = get_non_terminal_node_names(&rules, prefix);
            node_ids_by_name.extend(non_terminal_node_names);
            write_compressed_lines(&file, &rules, &encoded_paths, node_ids_by_name);
            log::info!("{} rules were written", rules.len());
        }
        Commands::Decompress { file, use_p_lines } => {
            log::info!("Decompressing file {}", file);
            let node_ids_by_name = parse_node_ids(&file, true);
            log::info!("Parsing lines");
            let (rules, paths) = parse_compressed_lines(&file, &node_ids_by_name);
            log::info!("Decompressing paths");
            let decompressed_paths = decompress(rules, paths);
            log::info!("Writing decompressed paths");
            write_decompressed_lines(&file, decompressed_paths, node_ids_by_name, use_p_lines);
        }
    }
}
