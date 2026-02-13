use crate::{
    helpers::{utils::NodeId, DeterministicHashMap, PathSegment, ReverseNodeRegistry},
    Rule,
};

pub fn print_grammar_simple(
    rules: &DeterministicHashMap<NodeId, Vec<NodeId>>,
    node_registry: &ReverseNodeRegistry,
) {
    for (meta_node, rule) in rules {
        let rule_text: String = rule
            .iter()
            .map(|n| node_registry.get_directed_name(*n))
            .collect();
        println!("Q\t{}\t{}", node_registry.get_name(meta_node.0), rule_text);
    }
}

pub fn print_grammar(grammar: &[Rule], node_registry: &ReverseNodeRegistry, print_addresses: bool) {
    for rule in grammar {
        if !print_addresses {
            println!(
                "Q\t{}\t{}{}",
                node_registry.get_name(rule.0 .0),
                node_registry.get_directed_name(rule.1),
                node_registry.get_directed_name(rule.2)
            );
        } else {
            println!(
                "Q\t{}\t{}{}\t{:?}",
                node_registry.get_name(rule.0 .0),
                node_registry.get_directed_name(rule.1),
                node_registry.get_directed_name(rule.2),
                rule.3
            );
        }
    }
}

pub fn print_walks(
    named_walks: &[(PathSegment, Vec<NodeId>)],
    node_registry: &ReverseNodeRegistry,
) {
    for (haplotype_name, walk) in named_walks {
        print!("W\t{}\t", haplotype_name);
        for node in walk {
            print!("{}", node_registry.get_directed_name(*node));
        }
        println!();
    }
}
