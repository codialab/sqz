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
    reverse_rule, Rule, Rules,
};
use itertools::Itertools;
use std::collections::{HashMap, HashSet};
use std::io::BufRead;

pub fn encode_paths(
    file: &str,
    rules: &Rules,
    offset: NodeId,
    node_ids_by_name: &HashMap<Vec<u8>, RawNodeId>,
) -> HashMap<String, Vec<NodeId>> {
    let mut result: HashMap<String, Vec<NodeId>> = HashMap::new();
    let mut rules_by_start_digram: HashMap<(NodeId, NodeId), Vec<&Rule>> = HashMap::new();
    for (_idx, rule) in rules.iter() {
        rules_by_start_digram
            .entry(canonize(rule.right[0], rule.right[1]))
            .and_modify(|e| e.push(rule))
            .or_insert(vec![rule]);
    }
    let mut rules_by_end_digram: HashMap<(NodeId, NodeId), Vec<&Rule>> = HashMap::new();
    for (_idx, rule) in rules.iter() {
        rules_by_end_digram
            .entry(canonize(
                rule.right[rule.right.len() - 2],
                rule.right[rule.right.len() - 1],
            ))
            .and_modify(|e| e.push(rule))
            .or_insert(vec![rule]);
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
            let encoded_path = encode_path(
                nodes,
                &rules_by_start_digram,
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

fn print_multiplicities(nodes: &Vec<NodeId>) {
    let mut nodes_visited_prev: HashSet<NodeId> = HashSet::new();
    let mut prev_counter = 0;
    let mut nodes_visited_curr: HashSet<NodeId> = HashSet::new();
    let mut curr_counter = 0;

    nodes.iter().tuple_windows().for_each(|(prev, curr)| {
        if nodes_visited_prev.contains(&prev.get_forward()) {
            prev_counter += 1;
            nodes_visited_prev.clear();
        }
        nodes_visited_prev.insert(prev.get_forward());
        if nodes_visited_curr.contains(&curr.get_forward()) {
            // log::error!("= Resetting outside curr =");
            curr_counter += 1;
            nodes_visited_curr.clear();
        }
        nodes_visited_curr.insert(curr.get_forward());
        print!("{}-({}|{})-", prev, prev_counter, curr_counter);
    });
    println!("{}", nodes.last().unwrap());
}

pub fn encode_path(
    nodes: Vec<NodeId>,
    rules_by_start: &HashMap<(NodeId, NodeId), Vec<&Rule>>,
    rules_by_end: &HashMap<(NodeId, NodeId), Vec<&Rule>>,
    statistics: &mut HashMap<NodeId, usize>,
    path_id: u64,
    offset: NodeId,
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
        // log::error!("Working on node {}, idx {}", current_node, idx);
        match stack.pop() {
            None => {
                stack.push(current_node);
                nodes_visited_prev.insert(current_node.get_forward());
            }
            Some(prev_node) => {
                if nodes_visited_curr.contains(&current_node.get_forward()) {
                    // log::error!("= Resetting outside curr =");
                    curr_counter += 1;
                    nodes_visited_curr.clear();
                }
                nodes_visited_curr.insert(current_node.get_forward());

                let (first_counter, second_counter) = (prev_counter, curr_counter);
                let (has_used_rule, _new_mult) = apply_rule(
                    &mut idx,
                    (prev_node, current_node),
                    (first_counter, second_counter),
                    &mut stack,
                    &mut multiplicity_stack,
                    &nodes,
                    statistics,
                    rules_by_start,
                    rules_by_end,
                    path_id,
                    &mut nodes_visited_prev,
                    &mut nodes_visited_curr,
                    &mut prev_counter,
                    &mut curr_counter,
                    offset,
                );

                if !has_used_rule {
                    // log::error!("Has not used rule for {} {}", prev_node, current_node);
                    stack.push(prev_node);
                    stack.push(current_node);
                    multiplicity_stack.push((first_counter, second_counter));
                } else if stack.len() >= 2 {
                    let mut change = true;
                    let mut counter = 0;
                    while change && stack.len() >= 2 {
                        let second = stack.pop().expect("stack has at least one node");
                        let first = stack.pop().expect("stack has at least two nodes");
                        let mult = multiplicity_stack
                            .pop()
                            .expect("multiplicity stack has at least one entry");
                        // Handle case beta := epsilon, but only if there is still sequence after it
                        let mult = if second > offset
                            && !second.is_forward()
                            && idx < nodes.len() - 1
                        {
                            log::error!("Special case of looking to the future (non-canonical), @ {}, orig: {:?}, new: {:?}", counter, mult, (mult.0, second_counter));
                            (mult.0, second_counter)
                        } else {
                            mult
                        };
                        log::error!("@ {}: using mult {:?}", counter, mult);
                        let result = apply_rule(
                            &mut idx,
                            (first, second),
                            mult,
                            &mut stack,
                            &mut multiplicity_stack,
                            &nodes,
                            statistics,
                            rules_by_start,
                            rules_by_end,
                            path_id,
                            &mut nodes_visited_prev,
                            &mut nodes_visited_curr,
                            &mut prev_counter,
                            &mut curr_counter,
                            offset,
                        );
                        change = result.0;
                        if !change {
                            stack.push(first);
                            stack.push(second);
                            multiplicity_stack.push(mult);
                        }
                        counter += 1;
                    }
                }

                // Do this last, to set prev only at the end and only once
                if nodes_visited_prev.contains(&current_node.get_forward()) {
                    // log::error!("= Resetting outside prev =");
                    prev_counter += 1;
                    nodes_visited_prev.clear()
                }
                nodes_visited_prev.insert(current_node.get_forward());
            }
        }
        // log::error!("stack: {:?}, idx: {}", stack, idx);
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
    rules_by_start: &HashMap<(NodeId, NodeId), Vec<&Rule>>,
    rules_by_end: &HashMap<(NodeId, NodeId), Vec<&Rule>>,
    path_id: u64,
    nodes_visited_prev: &mut HashSet<NodeId>,
    nodes_visited_curr: &mut HashSet<NodeId>,
    prev_counter: &mut usize,
    curr_counter: &mut usize,
    offset: NodeId,
) -> (bool, (usize, usize)) {
    // Handle case where we epsilon := beta
    let multiplicity =
        if digram.0 >= offset && digram.0.is_forward() && !multiplicity_stack.is_empty() {
            log::error!("Special case of looking to the past (canonical)");
            (
                multiplicity_stack[multiplicity_stack.len() - 1].1,
                multiplicity.1,
            )
        } else {
            multiplicity
        };
    let canonized_digram = canonize(digram.0, digram.1);
    let canonized_multiplicity = if is_edge_flipped(digram.0, digram.1) {
        transpose(multiplicity)
    } else {
        multiplicity
    };
    log::error!("Canonized: {:?}", canonized_digram);
    let mut has_used_rule = false;
    for (is_end, rule) in rules_by_start
        .get(&canonized_digram)
        .into_iter()
        .flatten()
        .map(|r| (false, *r))
        .chain(
            rules_by_end
                .get(&canonized_digram)
                .into_iter()
                .flatten()
                .map(|r| (true, *r)),
        )
    {
        log::error!(
            "\t-- Trying to work on rule: {}, with colors {:?}, mult: {:?} --",
            rule,
            rule.colors,
            multiplicity
        );
        if rule.colors.contains(
            path_id,
            canonized_multiplicity.0 as u32,
            canonized_multiplicity.1 as u32,
        ) {
            log::error!("\t== Working on rule: {} ==", rule);
            if (!is_end && digram.0 == rule.right[0])
                || (is_end && digram.0 == rule.right[rule.right.len() - 1].flip())
            {
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
                        "\tRemaining rule {:?} - {:?}, len first: {}, len second: {}, values {}, {}",
                        &rule.right[2..],
                        &nodes[*idx + 1..*idx + 1 + rule.right.len() - 2],
                        rule.right[2..].len(),
                        nodes[*idx + 1..*idx + 1 + rule.right.len() - 2].len(),
                        *idx, rule.right.len(),
                    );
                }
                if (rule.right.len() > 2
                    && !is_end
                    && nodes.len() > *idx + rule.right.len() - 2
                    && rule.right[2..] == nodes[*idx + 1..*idx + 1 + rule.right.len() - 2])
                    || (rule.right.len() > 2
                        && is_end
                        && nodes.len() - *idx > 1
                        && rule.right[..rule.right.len() - 2]
                            == nodes[*idx + 1..*idx + 1 + rule.right.len() - 2])
                    || (rule.right.len() == 2)
                {
                    // Forward match
                    if !is_end {
                        stack.push(rule.left.clone());
                    } else {
                        stack.push(rule.left.flip());
                    }
                    log::error!("\tInner forward:");
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
            } else if (!is_end && digram.0 == rule.right[1].flip())
                || (is_end && digram.0 == rule.right[rule.right.len() - 2])
            {
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
                        "Remaining rule {:?} - {:?}, is_end: {}, stack_len: {}, {:?} - {:?}",
                        &rule.right[2..],
                        reverse_rule(&stack[stack.len() - rule_remaining_len..].to_vec()),
                        is_end,
                        stack.len() >= rule_remaining_len,
                        &rule.right[..rule_remaining_len],
                        &stack[stack.len() - rule_remaining_len..]
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
                "\tNo color match: color: {:?}, {} | {} | {}",
                rule.colors,
                path_id,
                multiplicity.0,
                multiplicity.1
            );
        }
    }
    (has_used_rule, multiplicity)
}
