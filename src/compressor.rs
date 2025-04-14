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
    let mut counter_first_level = 0;
    let mut counter_all_levels = 0;
    let mut counter_total = 0;
    for (key, usage) in &statistics {
        if rules[key].colors.len() != *usage {
            counter_all_levels += 1;
            counter_total += (*usage as i64 - rules[key].colors.len() as i64).abs();
            let mut has_invalid_subrule = false;
            for subrule in &rules[key].right {
                if *subrule >= offset && rules[&subrule.get_forward()].colors.len() != statistics[&subrule.get_forward()] {
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
    panic!("Panicing");
    log::info!("Incorrect first level: {}, incorrect all levels: {}, incorrect total usages: {}", counter_first_level, counter_all_levels, counter_total);
    result
}

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
            .and_modify(|e| { panic!("Rule with same digram") })
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
            let encoded_path = encode_path2(
                nodes,
                &rules_by_digram,
                &mut statistics,
                path_id,
                offset,
            );
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
                if *subrule >= offset && rules[&subrule.get_forward()].colors.len() != statistics[&subrule.get_forward()] {
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
    log::info!("Incorrect first level: {}, incorrect all levels: {}, incorrect total usages: {}", counter_first_level, counter_all_levels, counter_total);
    result
}

pub fn encode_path2(
    nodes: Vec<NodeId>,
    rules_by_digram: &HashMap<(NodeId, NodeId), &Rule>,
    statistics: &mut HashMap<NodeId, usize>,
    path_id: u64,
    offset: NodeId,
) -> Vec<NodeId> {
    let mut stack: Vec<NodeId> = Vec::new();
    let mut multiplicity_stack: Vec<(usize, usize)> = Vec::new();

    let mut prev_counter = 0;
    let mut nodes_visited_curr: HashSet<NodeId> = HashSet::new();
    let mut curr_counter = 0;

    let mut node_idx = 0;

    let mut change_prev = None;
    log::info!("Encoding path {}", path_id);

    loop {
        // log::error!("Mult: {:?}, stack: {:?}",)
        if stack.len() >= 2 && is_rule_applicable(&(stack[stack.len() - 2], stack[stack.len() - 1]), path_id, multiplicity_stack[multiplicity_stack.len() - 1], rules_by_digram) {
            // Remove everything
            let second = stack.pop().unwrap();
            let first = stack.pop().unwrap();
            let mult = multiplicity_stack.pop().unwrap();

            // Apply rule
            let canonized = canonize(first, second);
            let left = rules_by_digram[&canonized].left;
            *statistics
                .get_mut(&left)
                .expect("statistics contains rule") += 1;
            if canonized.0 == first {
                stack.push(left);

                if !multiplicity_stack.is_empty() {
                    let old_mult = multiplicity_stack[multiplicity_stack.len() - 1];
                    log::error!("\tWriting forward: {}, mult_stack: {:?}, stack: {:?}", old_mult.1, multiplicity_stack, stack);
                    change_prev = Some(old_mult.1);
                }
            } else {
                stack.push(left.flip());

                if !multiplicity_stack.is_empty() && node_idx < nodes.len() - 1 {
                    let old_mult = multiplicity_stack[multiplicity_stack.len() - 1];
                    log::error!("\tWriting backward: {:?}", (old_mult.0, mult.1));
                    *multiplicity_stack.last_mut().expect("Mult stack has at least 1 element") = (old_mult.0, mult.1);
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
            if nodes_visited_curr.contains(&curr.get_forward()) {
                curr_counter += 1;
                nodes_visited_curr.clear();
            }
            nodes_visited_curr.insert(curr.get_forward());

            // Push to stack
            stack.push(curr);
            node_idx += 1;
            if stack.len() >= 2 {
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
    log::error!("Looking at {:?}, multiplicity {:?} (canon: {:?})", digram, multiplicity, canonized_multiplicity);
    if let Some(rule) = rules_by_digram.get(&canonized) {
        log::error!("\tMatching rule: {} -> {:?}, colors: {:?}", rule.left, rule.right, rule.colors);
        rule.colors.colors.contains(path_id, canonized_multiplicity.0 as u32, canonized_multiplicity.1 as u32)
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
                nodes_visited_curr.insert(current_node.get_forward());
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

                            log::error!("Special case of looking to the future (non-canonical), @ {}, orig: {:?}, new: {:?}, changed to {:?}", counter, mult, (mult.0, second_counter), mult);
                            (mult.0, second_counter)
                            // mult
                        } else {
                            mult
                        };
                        // log::error!("@ {}: using mult {:?}", counter, mult);
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

fn get_should_print_fully_generic(digram: &(NodeId, NodeId), printables: Vec<RawNodeId>) -> bool {
    for node in printables {
        if digram.0.get_forward() == NodeId::new(node - 1, 0) || digram.1.get_forward() == NodeId::new(node - 1, 0) {
            return true;
        }
    }
    false
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
    let should_print_nodes = vec![2619, 2637];
    // let should_print = get_should_print_fully_generic(&digram, should_print_nodes);
    let should_print = true;
    let canonized_digram = canonize(digram.0, digram.1);

    let multiplicity =
        if digram.0 >= offset && digram.0.is_forward() && !multiplicity_stack.is_empty() {
            if should_print {
                log::error!("Special case of looking to the past (canonical), orig: {:?}, new: {:?}", multiplicity, (multiplicity_stack[multiplicity_stack.len() - 1].1, multiplicity.1));
            }
            (
                multiplicity_stack[multiplicity_stack.len() - 1].1,
                multiplicity.1,
            )
        } else {
            multiplicity
        };
    let canonized_multiplicity = if is_edge_flipped(digram.0, digram.1) {
        transpose(multiplicity)
    } else {
        multiplicity
    };
    if should_print {
        log::error!("Canonized: {:?}, digram: {:?}, mult: {:?}, prev_c: {}, curr_c: {}", canonized_digram, digram, multiplicity, prev_counter, curr_counter);
        log::error!("Stack: {}, {:?}", stack.len(), stack.iter().rev().collect::<Vec<_>>());
        log::error!("Multiplicity stack: {}, {:?}", multiplicity_stack.len(), multiplicity_stack.iter().rev().collect::<Vec<_>>());
    }
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
        if should_print {
            log::error!(
                "\t-- Trying to work on rule: {}, with colors {:?}, mult: {:?}, path: {} --",
                rule,
                rule.colors,
                multiplicity,
                path_id,
            );
        }
        if rule.colors.colors.contains(
            path_id,
            canonized_multiplicity.0 as u32,
            canonized_multiplicity.1 as u32,
        ) {
            if should_print {
                // log::error!("\t== Working on rule: {} ==", rule);
            }
            if (!is_end && digram.0 == rule.right[0])
                || (is_end && digram.0 == rule.right[rule.right.len() - 1].flip())
            {
                // Exact match, continue forward
                if should_print {
                    // log::error!(
                    //     "{:?} - forward match | nodes_len {} | idx {} | rule_len {}",
                    //     digram,
                    //     nodes.len(),
                    //     idx,
                    //     rule.right.len()
                    // );
                }
                if rule.right.len() > 2 && nodes.len() > *idx + rule.right.len() - 2 {
                    if should_print {
                    // log::error!(
                    //     "\tRemaining rule {:?} - {:?}, len first: {}, len second: {}, values {}, {}",
                    //     &rule.right[2..],
                    //     &nodes[*idx + 1..*idx + 1 + rule.right.len() - 2],
                    //     rule.right[2..].len(),
                    //     nodes[*idx + 1..*idx + 1 + rule.right.len() - 2].len(),
                    //     *idx, rule.right.len(),
                    // );
                    }
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
                    if should_print {
                    // log::error!("\tInner forward:");
                    }
                    *statistics
                        .get_mut(&rule.left)
                        .expect("statistics contains rule") += 1;
                    has_used_rule = true;

                    let new_idx = *idx + rule.right.len() - 2;
                    for i in *idx + 1..=new_idx {
                        if nodes_visited_curr.contains(&nodes[i].get_forward()) {
                            // log::error!("= Resetting inside curr =");
                            *curr_counter += 1;
                            nodes_visited_curr.clear();
                        }
                        nodes_visited_curr.insert(nodes[i].get_forward());
                        if nodes_visited_prev.contains(&nodes[i].get_forward()) {
                            // log::error!("= Resetting inside prev =");
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
                if should_print {
                    // log::error!(
                    //     "{:?} - backward match, is_end: {},  rule_remaining_len {}, stack {} -> {} {}",
                    //     digram,
                    //     is_end,
                    //     rule_remaining_len,
                    //     stack.len(),
                    //     rule_remaining_len > 0,
                    //     stack.len() >= 1 + rule_remaining_len
                    // );
                }
                if rule_remaining_len > 0 && stack.len() >= rule_remaining_len {
                    if should_print {
                        // log::error!(
                        //     "Remaining rule {:?} - {:?}, is_end: {}, stack_len: {}, {:?} - {:?}",
                        //     &rule.right[2..],
                        //     reverse_rule(&stack[stack.len() - rule_remaining_len..].to_vec()),
                        //     is_end,
                        //     stack.len() >= rule_remaining_len,
                        //     &rule.right[..rule_remaining_len],
                        //     &stack[stack.len() - rule_remaining_len..]
                        // );
                    }
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
                    if should_print {
                        // log::error!("\tInner backward");
                    }
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
            if should_print {
                // log::error!(
                //     "\tNo color match: color: {:?}, {} | {} | {}",
                //     rule.colors,
                //     path_id,
                //     multiplicity.0,
                //     multiplicity.1
                // );
            }
        }
    }
    (has_used_rule, multiplicity)
}
