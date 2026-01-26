use itertools::Itertools;

use crate::{
    helpers::{
        utils::{NodeId, Orientation},
        PathSegment, ReverseNodeRegistry,
    },
    parser::Grammar,
};

pub fn decode_and_print_walks(
    haplotypes: Vec<(PathSegment, Vec<NodeId>)>,
    grammar: &Grammar,
    node_reg: &ReverseNodeRegistry,
    should_use_p_lines: bool,
) {
    for (name, walk) in haplotypes {
        let walk = decode_walk(walk, grammar);
        print_walk(name, walk, node_reg, should_use_p_lines);
    }
}

fn print_walk(
    name: PathSegment,
    walk: Vec<NodeId>,
    node_reg: &ReverseNodeRegistry,
    should_use_p_lines: bool,
) {
    if should_use_p_lines {
        print!("P\t{}\t", name.to_path_string());
    } else {
        print!("W\t{}\t", name.to_walk_string());
    }

    // Print all nodes except last
    for node in &walk[..walk.len() - 1] {
        if should_use_p_lines {
            let undirected = node.get_undirected();
            let name = node_reg.get_name(undirected);
            let direction = match node.1 {
                Orientation::Forward => '+',
                Orientation::Backward => '-',
            };
            print!("{}{},", name, direction);
        } else {
            let name = node_reg.get_directed_name(*node);
            print!("{}", name);
        }
    }

    // Print last node (necessary due to P-lines)
    let node = walk[walk.len() - 1];
    if should_use_p_lines {
        let undirected = node.get_undirected();
        let name = node_reg.get_name(undirected);
        let direction = match node.1 {
            Orientation::Forward => '+',
            Orientation::Backward => '-',
        };
        print!("{}{}", name, direction);
    } else {
        let name = node_reg.get_directed_name(node);
        print!("{}", name);
    }
    println!();
}

pub fn decode_walk(walk: Vec<NodeId>, grammar: &Grammar) -> Vec<NodeId> {
    let mut decoded: Vec<NodeId> = Vec::new();
    for node in walk {
        decoded.extend(decode_node(node, grammar));
    }
    decoded
}

fn flip_seq(seq: Vec<NodeId>) -> Vec<NodeId> {
    seq.into_iter().rev().map(|node| node.flip()).collect_vec()
}

fn decode_node(node: NodeId, grammar: &Grammar) -> Vec<NodeId> {
    if !node.is_meta_node() {
        return vec![node];
    }

    let rule = grammar
        .get(&node.get_undirected())
        .expect("All meta nodes have rules");
    let rule = vec![rule.0, rule.1];
    if node.is_forward() {
        rule.into_iter()
            .map(|n| decode_node(n, grammar))
            .flatten()
            .collect_vec()
    } else {
        flip_seq(rule)
            .into_iter()
            .map(|n| decode_node(n, grammar))
            .flatten()
            .collect_vec()
    }
}
