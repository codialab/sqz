use crate::{
    canonize,
    color_set::transpose,
    is_edge_flipped,
    node_id::NodeId,
    node_id::RawNodeId,
    parser::{
        bufreader_from_compressed_gfa, get_nodes_path, get_nodes_walk, parse_path_identifier,
        parse_walk_identifier,
    },
    Rule, Rules,
};
use itertools::Itertools;
use std::collections::{HashMap, HashSet};
use std::io::BufRead;

pub fn encode_paths2(
    file: &str,
    rules: &Rules,
    offset: NodeId,
    node_ids_by_name: &HashMap<Vec<u8>, RawNodeId>,
) -> HashMap<String, Vec<NodeId>> {
    let mut result: HashMap<String, Vec<NodeId>> = HashMap::new();
    let mut rules_by_digram: HashMap<(NodeId, NodeId), &Rule> = HashMap::new();
    for (_idx, rule) in rules.iter() {
        rules_by_digram
            .entry(canonize(rule.right[0], rule.right[1]))
            .and_modify(|_| panic!("Rule with same digram"))
            .or_insert(rule);
    }
    let mut data = bufreader_from_compressed_gfa(file);
    let mut statistics: HashMap<NodeId, usize> =
        rules.iter().map(|(k, _v)| (k.clone(), 0)).collect();
    let mut path_id = 0;

    // log::error!("RBD: {:?}", rules_by_digram);

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

            print_multiplicities(&nodes);
            let encoded_path =
                encode_path2(nodes, &rules_by_digram, &mut statistics, path_id, offset);
            result.insert(path_seg.to_string(), encoded_path);

            path_id += 1;
        }
        buf.clear();
    }
    let mut counter_first_level = 0;
    let mut counter_all_levels = 0;
    let mut counter_total = 0;
    for (key, usage) in &statistics {
        if rules[key].colors.len() != *usage {
            counter_all_levels += 1;
            counter_total += (*usage as i64 - rules[key].colors.len() as i64).abs();
            let mut has_invalid_subrule = false;
            for subrule in &rules[key].right {
                if *subrule >= offset
                    && rules[&subrule.get_forward()].colors.len()
                        != statistics[&subrule.get_forward()]
                {
                    has_invalid_subrule = true;
                }
            }
            if !has_invalid_subrule {
                counter_first_level += 1;
                log::warn!(
                    "Used {:?} -> {:?} exactly {} times vs supposed {} (offset: {})",
                    key,
                    rules[key].right,
                    usage,
                    rules[key].colors.len(),
                    offset
                );
            }
            // break;
        }
    }
    log::warn!(
        "Incorrect first level: {}, incorrect all levels: {}, incorrect total usages: {}",
        counter_first_level,
        counter_all_levels,
        counter_total
    );
    result
}

pub fn encode_path2(
    nodes: Vec<NodeId>,
    rules_by_digram: &HashMap<(NodeId, NodeId), &Rule>,
    statistics: &mut HashMap<NodeId, usize>,
    path_id: u64,
    _offset: NodeId,
) -> Vec<NodeId> {
    let mut stack: Vec<NodeId> = Vec::new();
    let mut multiplicity_stack: Vec<(usize, usize)> = Vec::new();

    let mut prev_counter;
    let mut nodes_visited_curr: HashSet<NodeId> = HashSet::new();
    let mut curr_counter = 0;

    let mut node_idx = 0;

    let mut change_prev = None;
    log::info!("Encoding path {}", path_id);

    loop {
        // log::error!("Mult: {:?}, stack: {:?}", multiplicity_stack, stack);
        if stack.len() >= 2
            && is_rule_applicable(
                &(stack[stack.len() - 2], stack[stack.len() - 1]),
                path_id,
                multiplicity_stack[multiplicity_stack.len() - 1],
                rules_by_digram,
            )
        {
            // Remove everything
            let second = stack.pop().unwrap();
            let first = stack.pop().unwrap();
            let mult = multiplicity_stack.pop().unwrap();

            // Apply rule
            let canonized = canonize(first, second);
            let left = rules_by_digram[&canonized].left;
            *statistics.get_mut(&left).expect("statistics contains rule") += 1;
            if canonized.0 == first {
                stack.push(left);

                if !multiplicity_stack.is_empty() {
                    let old_mult = multiplicity_stack[multiplicity_stack.len() - 1];
                    log::debug!("\tWriting forward: {}, mult_stack: {:?}, stack: {:?}", old_mult.1, multiplicity_stack, stack);
                    change_prev = Some(old_mult.1);
                } else {
                    log::debug!("Wanting to write forward, but cannot");
                }
            } else {
                stack.push(left.flip());

                if !multiplicity_stack.is_empty() {
                    let old_mult = multiplicity_stack[multiplicity_stack.len() - 1];
                    log::debug!("\tWriting backward: {:?}", (old_mult.0, mult.1));
                    *multiplicity_stack
                        .last_mut()
                        .expect("Mult stack has at least 1 element") = (old_mult.0, mult.1);
                } else if !multiplicity_stack.is_empty() {
                    let old_mult = multiplicity_stack[multiplicity_stack.len() - 1];
                    log::debug!("Writing backward over the end, er: {}, mult: {:?}, ", curr_counter, mult);
                    *multiplicity_stack
                        .last_mut()
                        .expect("Mult stack has at least 1 element") = (old_mult.0, curr_counter);
                }
            }
        } else if node_idx < nodes.len() {
            let curr = nodes[node_idx];
            // Update counters
            if change_prev.is_some() {
                prev_counter = change_prev.unwrap();
                change_prev = None;
            } else {
                prev_counter = curr_counter;
            }
            if !stack.is_empty() {
                curr_counter += 1;
            }
            // if nodes_visited_curr.contains(&curr.get_forward()) {
            //     curr_counter += 1;
            //     nodes_visited_curr.clear();
            // }
            // nodes_visited_curr.insert(curr.get_forward());

            // Push to stack
            stack.push(curr);
            node_idx += 1;
            if stack.len() >= 1 {
                multiplicity_stack.push((prev_counter, curr_counter))
            }
        } else {
            break;
        }
    }
    stack
}

fn is_rule_applicable(
    digram: &(NodeId, NodeId),
    path_id: u64,
    multiplicity: (usize, usize),
    rules_by_digram: &HashMap<(NodeId, NodeId), &Rule>,
) -> bool {
    let canonized = canonize(digram.0, digram.1);
    let canonized_multiplicity = if is_edge_flipped(digram.0, digram.1) {
        transpose(multiplicity)
    } else {
        multiplicity
    };
    log::debug!(
        "Looking at {:?}, multiplicity {:?} (canon: {:?})",
        digram,
        multiplicity,
        canonized_multiplicity
    );
    if let Some(rule) = rules_by_digram.get(&canonized) {
        log::debug!(
            "\tMatching rule: {} -> {:?}, colors: {:?}",
            rule.left,
            rule.right,
            rule.colors
        );
        log::debug!(
            "\tMatched: {}",
            rule.colors.colors.contains(
                path_id,
                canonized_multiplicity.0 as u32,
                canonized_multiplicity.1 as u32
            )
        );
        rule.colors.colors.contains(
            path_id,
            canonized_multiplicity.0 as u32,
            canonized_multiplicity.1 as u32,
        )
    } else {
        false
    }
}

fn print_multiplicities(nodes: &Vec<NodeId>) {
    let mut nodes_visited_prev: HashSet<NodeId> = HashSet::new();
    let mut prev_counter = 0;
    let mut nodes_visited_curr: HashSet<NodeId> = HashSet::new();
    let mut curr_counter = 0;

    nodes_visited_curr.insert(nodes[0].get_forward());
    nodes.iter().tuple_windows().for_each(|(prev, curr)| {
        prev_counter = curr_counter;
        curr_counter += 1;
        // if nodes_visited_prev.contains(&prev.get_forward()) {
        //     prev_counter += 1;
        //     nodes_visited_prev.clear();
        // }
        // nodes_visited_prev.insert(prev.get_forward());
        // if nodes_visited_curr.contains(&curr.get_forward()) {
        //     // log::error!("= Resetting outside curr =");
        //     curr_counter += 1;
        //     nodes_visited_curr.clear();
        // }
        nodes_visited_curr.insert(curr.get_forward());
        print!("{}-({}|{})-", prev, prev_counter, curr_counter);
    });
    println!("{}", nodes.last().unwrap());
}
