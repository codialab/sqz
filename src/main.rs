use anyhow::Result;
use clap::{Parser, Subcommand};
use helpers::{
    utils::Address, utils::AddressNumber, utils::CanonicalDigram, utils::Digram, digram_occurrences::DigramOccurrences,
    utils::LocalizedDigram, utils::NodeId, NodeRegistry, Occurrence,
};
use itertools::Itertools;
use parser::parse_file_to_haplotypes;
use std::{
    collections::{HashMap, HashSet},
    hash::{BuildHasherDefault, DefaultHasher},
    path::PathBuf,
};

use crate::helpers::{PathSegment, ReverseNodeRegistry};

mod helpers;
mod parser;

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
        file: PathBuf,

        /// Prefix that is used for non-terminal node identifiers
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

type Rule = (
    NodeId,
    NodeId,
    NodeId,
    HashSet<Occurrence, BuildHasherDefault<DefaultHasher>>,
);

fn get_canonized(
    u: NodeId,
    v: NodeId,
    haplotype: usize,
    a_0: AddressNumber,
    a_1: AddressNumber,
) -> (CanonicalDigram, Occurrence) {
    let digram = Digram::new(u, v);
    let address = Address::from_address_number(a_0, a_1);
    let (digram, address) = LocalizedDigram::new(digram, address).split_to_canonical();
    let occurrence = Occurrence::new(haplotype, address);
    (digram, occurrence)
}

fn get_haplotype_walks(d: &DigramOccurrences, rules: &Vec<Rule>, singleton_haplotypes: &HashMap<usize, NodeId>, number_of_paths: usize) -> Vec<Vec<NodeId>> {
    let mut walks: Vec<Vec<(Digram, Address)>> = vec![Vec::new(); number_of_paths];

    // Find all singleton digrams
    for (digram, occurrences) in d {
        for occurrence in occurrences {
            let walk_id = occurrence.get_haplotype();
            let address = occurrence.get_address();
            if address.is_forward() {
                walks[walk_id].push((digram.clone().into(), address));
            } else {
                walks[walk_id].push((Into::<Digram>::into(digram.clone()).flip(), address.flip()));
            }
        }
    }

    // Sort all singleton digrams
    let mut walks = walks.into_iter().map(|mut walk| {
        walk.sort_by(|a, b| a.1.cmp(&b.1));
        if walk.len() > 0 {
            let first_value = vec![walk[0].0.0];
            first_value.into_iter().chain(walk.into_iter().map(|(digram, _address)| digram.1)).collect_vec()
        } else {
            Vec::new()
        }
    }).collect_vec();

    // Insert all walks that only consist of a single node
    let single_node_walks = walks.iter().enumerate().filter_map(|(idx, walk)| if walk.is_empty() { Some(idx) } else { None }).collect_vec();
    log::info!("SNW: {:?}", single_node_walks);
    for single_node_walk in single_node_walks {
        if singleton_haplotypes.contains_key(&single_node_walk) {
            walks[single_node_walk].push(singleton_haplotypes[&single_node_walk]);
        } else {
            if let Some(meta_node) = get_single_meta_node_walk(single_node_walk, rules) {
                walks[single_node_walk].push(meta_node);
            } else {
                log::error!("Was not able to assign a sequence to walk No. {}", single_node_walk);
            }
        }
    }
    walks
}

fn get_single_meta_node_walk(single_node_walk: usize, rules: &Vec<Rule>) -> Option<NodeId> {
    for rule in rules.iter().rev() {
        for occurrence in &rule.3 {
            if occurrence.get_haplotype() == single_node_walk {
                if occurrence.get_address().is_forward() {
                    return Some(rule.0);
                } else {
                    return Some(rule.0.flip());
                }
            }
        }
    }
    None
}

fn build_rule(
    uv: CanonicalDigram,
    d: &mut DigramOccurrences,
    node_registry: &mut NodeRegistry,
) -> Rule {
    let q = node_registry.get_new_meta_node();
    let mut d_q = d.remove_digram(&uv);
    for occurrence in &d_q {
        let Some((left_neighbor, left_address_number)) = d.get_left_neighbor(&uv, occurrence)
        else {
            continue;
        };
        let (old_digram, old_occurrence) = get_canonized(
            left_neighbor,
            uv.get_u(),
            occurrence.get_haplotype(),
            left_address_number,
            occurrence.get_address().get_first(),
        );
        d.delete_occurrence(&old_digram, &old_occurrence);
        let (new_digram, new_occurrence) = get_canonized(
            left_neighbor,
            q,
            occurrence.get_haplotype(),
            left_address_number,
            occurrence.get_address().get_first(),
        );
        d.add_occurrence(&new_digram, new_occurrence);
    }
    if uv.get_u() == uv.get_v() {
        let (new_d_q, d_qq, d_qv) = Occurrence::split_self_loops(d_q);
        println!("Self loop:");
        println!("\td_qq: {:?}", d_qq);
        println!("\td_q: {:?}", new_d_q);
        println!("\td_qv: {:?}", d_qv);
        d_q = new_d_q.into_iter().collect();
        d.add_digram(&Digram::new(q, q).into(), d_qq.into_iter().collect());
        d.add_digram(
            &Digram::new(q, uv.get_v()).into(),
            d_qv.into_iter().collect(),
        );
    }
    // } else if d.contains(&Digram::new(uv.get_v(), uv.get_u()).into()) {
    //     d.remove_digram(&Digram::new(uv.get_v(), uv.get_u()).into())
    //         .into_iter()
    //         .collect()
    // } else {
    //     Vec::new()
    // };
    // if !d_vu.is_empty() {
    //     let (d_qq, d_vq, d_qu, d_vu) = Occurrence::split_pre_self_loop(d_vu, &d_q);
    // }
    for occurrence in &d_q {
        let Some((right_neighbor, right_address_number)) = d.get_right_neighbor(&uv, occurrence)
        else {
            continue;
        };
        let (old_digram, old_occurrence) = get_canonized(
            uv.get_v(),
            right_neighbor,
            occurrence.get_haplotype(),
            occurrence.get_address().get_second(),
            right_address_number,
        );
        d.delete_occurrence(&old_digram, &old_occurrence);
        let (new_digram, new_occurrence) = get_canonized(
            q,
            right_neighbor,
            occurrence.get_haplotype(),
            occurrence.get_address().get_first(),
            right_address_number,
        );
        d.add_occurrence(&new_digram, new_occurrence);
    }
    (q, uv.get_u(), uv.get_v(), d_q)
}

fn build_grammar(d: &mut DigramOccurrences, node_registry: &mut NodeRegistry) -> Vec<Rule> {
    let mut rules: Vec<Rule> = Vec::new();
    while let Some(uv) = d.get_most_frequent(2) {
        rules.push(build_rule(uv, d, node_registry));
    }
    rules
}

fn print_grammar(grammar: &Vec<Rule>, node_registry: &ReverseNodeRegistry, print_addresses: bool) {
    for rule in grammar {
        if !print_addresses {
            println!("Q\t{}\t{}{}", node_registry.get_name(rule.0.0), node_registry.get_directed_name(rule.1), node_registry.get_directed_name(rule.2));
        } else {
            println!("Q\t{}\t{}{}\t{:?}", node_registry.get_name(rule.0.0), node_registry.get_directed_name(rule.1), node_registry.get_directed_name(rule.2), rule.3);
        }
    } 
}

fn print_walks(walks: &Vec<Vec<NodeId>>, node_registry: &ReverseNodeRegistry, haplotype_names: &Vec<PathSegment>) {
    for (walk, haplotype_name) in walks.iter().zip(haplotype_names.iter()) {
        print!("W\t{}\t", haplotype_name);
        for node in walk {
            print!("{}", node_registry.get_directed_name(*node));
        }
        println!("");
    }
}

fn main() -> Result<()> {
    env_logger::init();
    let args = Cli::parse();

    match args.command {
        Commands::Compress { file, prefix: _ } => {
            log::info!("Compressing file {:?}", file);
            let (haplotypes, mut node_registry, singleton_haplotypes) = parse_file_to_haplotypes(&file, true)?;
            let (haplotype_names, haplotypes): (Vec<PathSegment>, Vec<Vec<LocalizedDigram>>) = haplotypes.into_iter().unzip(); 
            let number_of_paths = haplotypes.len();
            log::info!("Parsed {} haplotypes", number_of_paths);
            let mut d = DigramOccurrences::from(haplotypes);
            let grammar = build_grammar(&mut d, &mut node_registry);
            let rev_reg: ReverseNodeRegistry = node_registry.into();
            for (digram, occurrences) in &d {
                log::info!("| {}{} -> {:?}", rev_reg.get_directed_name(digram.0), rev_reg.get_directed_name(digram.1), occurrences);
            }
            print_grammar(&grammar, &rev_reg, true);
            let haplotype_walks = get_haplotype_walks(&d, &grammar, &singleton_haplotypes, number_of_paths);
            print_walks(&haplotype_walks, &rev_reg, &haplotype_names);
        }
        Commands::Decompress { file, use_p_lines: _ } => {
            log::info!("Decompressing file {}", file);
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use crate::parser::{get_haplotypes_from_walk_strings};

    use super::*;

    #[test]
    fn test_loop_free_rule_creation() {
        let haplotypes = vec![">0>1>2>3", ">4>1>2>3", ">5>1>2>6"];
        let mut node_registry = NodeRegistry::new();
        let haplotypes = get_haplotypes_from_walk_strings(haplotypes, &mut node_registry);
        let mut d = DigramOccurrences::from(haplotypes);
        let u = NodeId::new(
            node_registry.get_id("1".as_bytes()),
            helpers::utils::Orientation::Forward,
        );
        let v = NodeId::new(
            node_registry.get_id("2".as_bytes()),
            helpers::utils::Orientation::Forward,
        );
        let uv = Digram::new(u, v).into();
        assert!(!d.are_neighbors_equal());
        let rule = build_rule(uv, &mut d, &mut node_registry);
        assert!(!d.are_neighbors_equal());
        assert_eq!(rule.1, u);
        assert_eq!(rule.2, v);
        assert_eq!(rule.3.len(), 3);
    }

    #[test]
    fn test_loop_free_grammar_creation() {
        let haplotypes = vec![">0>1>2>3", ">4>1>2>3", ">5>1>2>6"];
        let mut node_registry = NodeRegistry::new();
        let haplotypes = get_haplotypes_from_walk_strings(haplotypes, &mut node_registry);
        println!("haplotypes: {:?}", haplotypes);
        let mut d = DigramOccurrences::from(haplotypes);
        println!("d: {:?}", d);
        let grammar = build_grammar(&mut d, &mut node_registry);
        assert_eq!(grammar.len(), 2);
    }

    #[test]
    fn test_self_loop_rule_creation() {
        let haplotypes = vec![">0>1>1>1>1>2"];
        let mut node_registry = NodeRegistry::new();
        let haplotypes = get_haplotypes_from_walk_strings(haplotypes, &mut node_registry);
        let mut d = DigramOccurrences::from(haplotypes);
        let u = NodeId::new(
            node_registry.get_id("1".as_bytes()),
            helpers::utils::Orientation::Forward,
        );
        let v = NodeId::new(
            node_registry.get_id("1".as_bytes()),
            helpers::utils::Orientation::Forward,
        );
        let uv = Digram::new(u, v).into();
        assert!(!d.are_neighbors_equal());
        let rule = build_rule(uv, &mut d, &mut node_registry);
        println!("Rule: {:?}", rule);
        assert!(!d.are_neighbors_equal());
        assert_eq!(rule.1, u);
        assert_eq!(rule.2, v);
        assert_eq!(rule.3.len(), 2);
    }

    #[test]
    fn test_double_self_loop_rule_creation() {
        let haplotypes = vec![">0>1>1>1>1>1>1>1>1>2"];
        let mut node_registry = NodeRegistry::new();
        let haplotypes = get_haplotypes_from_walk_strings(haplotypes, &mut node_registry);
        let mut d = DigramOccurrences::from(haplotypes);
        let u = NodeId::new(
            node_registry.get_id("1".as_bytes()),
            helpers::utils::Orientation::Forward,
        );
        let v = NodeId::new(
            node_registry.get_id("1".as_bytes()),
            helpers::utils::Orientation::Forward,
        );
        let uv = Digram::new(u, v).into();
        assert!(!d.are_neighbors_equal());
        let rule = build_rule(uv, &mut d, &mut node_registry);
        println!("Rule: {:?}", rule);
        assert!(!d.are_neighbors_equal());
        assert_eq!(rule.1, u);
        assert_eq!(rule.2, v);
        assert_eq!(rule.3.len(), 4);
        let self_loop_meta =  rule.0;
        let uv = Digram::new(self_loop_meta, self_loop_meta).into();
        let rule = build_rule(uv, &mut d, &mut node_registry);
        println!("Rule: {:?}", rule);
        assert!(!d.are_neighbors_equal());
        assert_eq!(rule.1, self_loop_meta);
        assert_eq!(rule.2, self_loop_meta);
        assert_eq!(rule.3.len(), 2);
    }

    #[test]
    fn test_double_self_loop_rule_creation_reverse() {
        let haplotypes = vec!["<0<1<1<1<1<1<1<1<1<2"];
        let mut node_registry = NodeRegistry::new();
        let haplotypes = get_haplotypes_from_walk_strings(haplotypes, &mut node_registry);
        let mut d = DigramOccurrences::from(haplotypes);
        let u = NodeId::new(
            node_registry.get_id("1".as_bytes()),
            helpers::utils::Orientation::Backward,
        );
        let v = NodeId::new(
            node_registry.get_id("1".as_bytes()),
            helpers::utils::Orientation::Backward,
        );
        let uv = Digram::new(u, v).into();
        assert!(!d.are_neighbors_equal());
        let rule = build_rule(uv, &mut d, &mut node_registry);
        println!("Rule: {:?}", rule);
        assert!(!d.are_neighbors_equal());
        assert_eq!(rule.1, u.flip());
        assert_eq!(rule.2, v.flip());
        assert_eq!(rule.3.len(), 4);
        let self_loop_meta =  rule.0;
        let uv = Digram::new(self_loop_meta, self_loop_meta).into();
        let rule = build_rule(uv, &mut d, &mut node_registry);
        println!("Rule: {:?}", rule);
        assert!(!d.are_neighbors_equal());
        assert_eq!(rule.1, self_loop_meta);
        assert_eq!(rule.2, self_loop_meta);
        assert_eq!(rule.3.len(), 2);
    }

    #[test]
    fn test_one_half_self_loop_rule_creation() {
        let haplotypes = vec![">0>1>1>1>1>1>1>2"];
        let mut node_registry = NodeRegistry::new();
        let haplotypes = get_haplotypes_from_walk_strings(haplotypes, &mut node_registry);
        let mut d = DigramOccurrences::from(haplotypes);
        let u = NodeId::new(
            node_registry.get_id("1".as_bytes()),
            helpers::utils::Orientation::Forward,
        );
        let v = NodeId::new(
            node_registry.get_id("1".as_bytes()),
            helpers::utils::Orientation::Forward,
        );
        let uv = Digram::new(u, v).into();
        assert!(!d.are_neighbors_equal());
        let rule = build_rule(uv, &mut d, &mut node_registry);
        println!("Rule: {:?}", rule);
        assert!(!d.are_neighbors_equal());
        assert_eq!(rule.1, u);
        assert_eq!(rule.2, v);
        assert_eq!(rule.3.len(), 3);
        let self_loop_meta =  rule.0;
        let uv = Digram::new(self_loop_meta, self_loop_meta).into();
        let rule = build_rule(uv, &mut d, &mut node_registry);
        println!("Rule: {:?}", rule);
        assert!(!d.are_neighbors_equal());
        assert_eq!(rule.1, self_loop_meta);
        assert_eq!(rule.2, self_loop_meta);
        assert_eq!(rule.3.len(), 1);
    }

    #[test]
    fn test_one_half_self_loop_rule_creation_reverse() {
        let haplotypes = vec!["<0<1<1<1<1<1<1<2"];
        let mut node_registry = NodeRegistry::new();
        let haplotypes = get_haplotypes_from_walk_strings(haplotypes, &mut node_registry);
        let mut d = DigramOccurrences::from(haplotypes);
        let u = NodeId::new(
            node_registry.get_id("1".as_bytes()),
            helpers::utils::Orientation::Forward,
        );
        let v = NodeId::new(
            node_registry.get_id("1".as_bytes()),
            helpers::utils::Orientation::Forward,
        );
        let uv = Digram::new(u, v).into();
        assert!(!d.are_neighbors_equal());
        let rule = build_rule(uv, &mut d, &mut node_registry);
        println!("Rule: {:?}", rule);
        assert!(!d.are_neighbors_equal());
        assert_eq!(rule.1, u);
        assert_eq!(rule.2, v);
        assert_eq!(rule.3.len(), 3);
        let self_loop_meta =  rule.0;
        let uv = Digram::new(self_loop_meta, self_loop_meta).into();
        let rule = build_rule(uv, &mut d, &mut node_registry);
        println!("Rule: {:?}", rule);
        assert!(!d.are_neighbors_equal());
        assert_eq!(rule.1, self_loop_meta);
        assert_eq!(rule.2, self_loop_meta);
        assert_eq!(rule.3.len(), 1);
    }

    #[test]
    fn test_pre_self_loop_case_2() {
        let haplotypes = vec![">0>1>2>1>2>3"];
        let mut node_registry = NodeRegistry::new();
        let haplotypes = get_haplotypes_from_walk_strings(haplotypes, &mut node_registry);
        let mut d = DigramOccurrences::from(haplotypes);
        let u = NodeId::new(
            node_registry.get_id("1".as_bytes()),
            helpers::utils::Orientation::Forward,
        );
        let v = NodeId::new(
            node_registry.get_id("2".as_bytes()),
            helpers::utils::Orientation::Forward,
        );
        let uv = Digram::new(u, v).into();
        assert!(!d.are_neighbors_equal());
        let rule = build_rule(uv, &mut d, &mut node_registry);
        println!("Rule: {:?}", rule);
        assert!(!d.are_neighbors_equal());
        assert_eq!(rule.1, u);
        assert_eq!(rule.2, v);
        assert_eq!(rule.3.len(), 2);
        let self_loop_meta =  rule.0;
        let uv = Digram::new(self_loop_meta, self_loop_meta).into();
        let rule = build_rule(uv, &mut d, &mut node_registry);
        println!("Rule: {:?}", rule);
        assert!(!d.are_neighbors_equal());
        assert_eq!(rule.1, self_loop_meta);
        assert_eq!(rule.2, self_loop_meta);
        assert_eq!(rule.3.len(), 1);
    }

    #[test]
    fn test_pre_self_loop_case_2_reverse() {
        let haplotypes = vec!["<0<1<2<1<2<3"];
        let mut node_registry = NodeRegistry::new();
        let haplotypes = get_haplotypes_from_walk_strings(haplotypes, &mut node_registry);
        let mut d = DigramOccurrences::from(haplotypes);
        let u = NodeId::new(
            node_registry.get_id("1".as_bytes()),
            helpers::utils::Orientation::Backward,
        );
        let v = NodeId::new(
            node_registry.get_id("2".as_bytes()),
            helpers::utils::Orientation::Backward,
        );
        let uv = Digram::new(u, v).into();
        assert!(!d.are_neighbors_equal());
        let rule = build_rule(uv, &mut d, &mut node_registry);
        println!("Rule: {:?}", rule);
        assert!(!d.are_neighbors_equal());
        assert_eq!(rule.1, u);
        assert_eq!(rule.2, v);
        assert_eq!(rule.3.len(), 2);
        let self_loop_meta =  rule.0;
        let uv = Digram::new(self_loop_meta, self_loop_meta).into();
        let rule = build_rule(uv, &mut d, &mut node_registry);
        println!("Rule: {:?}", rule);
        assert!(!d.are_neighbors_equal());
        assert_eq!(rule.1, self_loop_meta);
        assert_eq!(rule.2, self_loop_meta);
        assert_eq!(rule.3.len(), 1);
    }

    #[test]
    fn test_pre_self_loop_case_3() {
        let haplotypes = vec![">0>1>2>1>3"];
        let mut node_registry = NodeRegistry::new();
        let haplotypes = get_haplotypes_from_walk_strings(haplotypes, &mut node_registry);
        let mut d = DigramOccurrences::from(haplotypes);
        let u = NodeId::new(
            node_registry.get_id("1".as_bytes()),
            helpers::utils::Orientation::Forward,
        );
        let v = NodeId::new(
            node_registry.get_id("2".as_bytes()),
            helpers::utils::Orientation::Forward,
        );
        let uv = Digram::new(u, v).into();
        assert!(!d.are_neighbors_equal());
        let rule = build_rule(uv, &mut d, &mut node_registry);
        println!("Rule: {:?}", rule);
        assert!(!d.are_neighbors_equal());
        assert_eq!(rule.1, u);
        assert_eq!(rule.2, v);
        assert_eq!(rule.3.len(), 1);
        let self_loop_meta =  rule.0;
        let uv = Digram::new(self_loop_meta, u).into();
        let rule = build_rule(uv, &mut d, &mut node_registry);
        println!("Rule: {:?}", rule);
        assert!(!d.are_neighbors_equal());
        assert_eq!(rule.1, u.flip());
        assert_eq!(rule.2, self_loop_meta.flip());
        assert_eq!(rule.3.len(), 1);
    }

    #[test]
    fn test_pre_self_loop_case_3_reverse() {
        let haplotypes = vec!["<0<1<2<1<3"];
        let mut node_registry = NodeRegistry::new();
        let haplotypes = get_haplotypes_from_walk_strings(haplotypes, &mut node_registry);
        let mut d = DigramOccurrences::from(haplotypes);
        let u = NodeId::new(
            node_registry.get_id("1".as_bytes()),
            helpers::utils::Orientation::Backward,
        );
        let v = NodeId::new(
            node_registry.get_id("2".as_bytes()),
            helpers::utils::Orientation::Backward,
        );
        let uv = Digram::new(u, v).into();
        assert!(!d.are_neighbors_equal());
        let rule = build_rule(uv, &mut d, &mut node_registry);
        println!("Rule: {:?}", rule);
        assert!(!d.are_neighbors_equal());
        assert_eq!(rule.1, u);
        assert_eq!(rule.2, v);
        assert_eq!(rule.3.len(), 1);
        let self_loop_meta =  rule.0;
        let uv = Digram::new(self_loop_meta, u).into();
        let rule = build_rule(uv, &mut d, &mut node_registry);
        println!("Rule: {:?}", rule);
        assert!(!d.are_neighbors_equal());
        assert_eq!(rule.1, u.flip());
        assert_eq!(rule.2, self_loop_meta.flip());
        assert_eq!(rule.3.len(), 1);
    }

    #[test]
    fn test_pre_self_loop_case_4() {
        let haplotypes = vec![">0>1>0>2>1>2>3"];
        let mut node_registry = NodeRegistry::new();
        let haplotypes = get_haplotypes_from_walk_strings(haplotypes, &mut node_registry);
        let mut d = DigramOccurrences::from(haplotypes);
        let u = NodeId::new(
            node_registry.get_id("1".as_bytes()),
            helpers::utils::Orientation::Forward,
        );
        let v = NodeId::new(
            node_registry.get_id("2".as_bytes()),
            helpers::utils::Orientation::Forward,
        );
        let uv = Digram::new(u, v).into();
        assert!(!d.are_neighbors_equal());
        let rule = build_rule(uv, &mut d, &mut node_registry);
        println!("Rule: {:?}", rule);
        assert!(!d.are_neighbors_equal());
        assert_eq!(rule.1, u);
        assert_eq!(rule.2, v);
        assert_eq!(rule.3.len(), 1);
        let self_loop_meta =  rule.0;
        let uv = Digram::new(v, self_loop_meta).into();
        let rule = build_rule(uv, &mut d, &mut node_registry);
        println!("Rule: {:?}", rule);
        assert!(!d.are_neighbors_equal());
        assert_eq!(rule.1, v);
        assert_eq!(rule.2, self_loop_meta);
        assert_eq!(rule.3.len(), 1);
    }

    #[test]
    fn test_pre_self_loop_case_4_reverse() {
        let haplotypes = vec!["<0<1<0<2<1<2<3"];
        let mut node_registry = NodeRegistry::new();
        let haplotypes = get_haplotypes_from_walk_strings(haplotypes, &mut node_registry);
        let mut d = DigramOccurrences::from(haplotypes);
        let u = NodeId::new(
            node_registry.get_id("1".as_bytes()),
            helpers::utils::Orientation::Backward,
        );
        let v = NodeId::new(
            node_registry.get_id("2".as_bytes()),
            helpers::utils::Orientation::Backward,
        );
        let uv = Digram::new(u, v).into();
        assert!(!d.are_neighbors_equal());
        let rule = build_rule(uv, &mut d, &mut node_registry);
        println!("Rule: {:?}", rule);
        assert!(!d.are_neighbors_equal());
        assert_eq!(rule.1, u);
        assert_eq!(rule.2, v);
        assert_eq!(rule.3.len(), 1);
        let self_loop_meta =  rule.0;
        let uv = Digram::new(v, self_loop_meta).into();
        let rule = build_rule(uv, &mut d, &mut node_registry);
        println!("Rule: {:?}", rule);
        assert!(!d.are_neighbors_equal());
        assert_eq!(rule.1, v);
        assert_eq!(rule.2, self_loop_meta);
        assert_eq!(rule.3.len(), 1);
    }
}
