use crate::{
    helpers::{utils::NodeId, PathSegment, ReverseNodeRegistry},
    Rule,
};

pub fn print_grammar(
    grammar: &Vec<Rule>,
    node_registry: &ReverseNodeRegistry,
    print_addresses: bool,
) {
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
    walks: &Vec<Vec<NodeId>>,
    node_registry: &ReverseNodeRegistry,
    haplotype_names: &Vec<PathSegment>,
) {
    for (walk, haplotype_name) in walks.iter().zip(haplotype_names.iter()) {
        print!("W\t{}\t", haplotype_name);
        for node in walk {
            print!("{}", node_registry.get_directed_name(*node));
        }
        println!("");
    }
}
