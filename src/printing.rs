use crate::helpers::{
    utils::{NodeId, UndirectedNodeId},
    DeterministicHashMap, PathSegment, ReverseNodeRegistry,
};

pub fn print_grammar_simple(
    rules: &DeterministicHashMap<UndirectedNodeId, Vec<NodeId>>,
    node_registry: &ReverseNodeRegistry,
) {
    for (meta_node, rule) in rules {
        let rule_text: String = rule
            .iter()
            .map(|n| node_registry.get_directed_name(*n))
            .collect();
        println!("Q\t{}\t{}", node_registry.get_name(*meta_node), rule_text);
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
