use anyhow::Result;
use clap::{Parser, Subcommand};
use deepsize::DeepSizeOf;
use env_logger::Env;
use helpers::{digram_occurrences::DigramOccurrences, utils::LocalizedDigram};
use itertools::Itertools;
use std::{
    collections::{HashMap, HashSet, VecDeque},
    hash::BuildHasherDefault,
    path::PathBuf,
};

use crate::{
    compressing::{compress, compress_remaining_file},
    decoding::decode_and_print_walks,
    encoding::get_haplotype_walks,
    grammar_building::build_grammar,
    helpers::{
        utils::{CanonicalDigram, Digram, NodeId, UndirectedNodeId},
        DeterministicHashMap, NodeRegistry, PathSegment, ReverseNodeRegistry,
    },
    parser::{
        compare_file, parse_file_to_digrams, parse_file_to_digrams_ratio_based,
        parse_file_to_haplotypes_with_grammar, Grammar,
    },
    printing::{print_grammar_simple, print_walks},
};

mod aho_corasick;
mod compressing;
mod decoding;
mod encoding;
mod grammar_building;
mod helpers;
mod parser;
mod printing;

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Debug, Subcommand)]
enum Commands {
    /// Compresses a GFA file by building the grammar on the full
    /// number of paths
    #[command(arg_required_else_help = true)]
    CompressFull {
        /// Input GFA file
        #[arg(required = true)]
        file: PathBuf,

        /// Prefix that is used for non-terminal node identifiers
        #[arg(short, long, default_value = "Q")]
        prefix: String,

        /// Do not reduce number of rules based on utility (keeps all rules
        /// at a length of 2)
        #[arg(short = 'n', long)]
        no_reduction: bool,
    },

    /// Decompresses a GFA file (by default paths will be decompressed as W-lines)
    #[command(arg_required_else_help = true)]
    Decompress {
        /// Input GFA file
        #[arg(required = true)]
        file: PathBuf,

        #[arg(short = 'p', long)]
        use_p_lines: bool,
    },

    /// Compresses a GFA file by building the grammar from a subset
    /// of paths.
    #[command(arg_required_else_help = true)]
    CompressPartial {
        /// Input GFA file
        #[arg(required = true)]
        file: PathBuf,

        /// Prefix that is used for non-terminal node identifiers
        #[arg(short, long, default_value = "Q")]
        prefix: String,

        /// Ratio of paths used for the grammar building
        #[arg(short, long, default_value = "0.5")]
        ratio: f32,

        /// Do not reduce number of rules based on utility (keeps all rules
        /// at a length of 2)
        #[arg(short = 'n', long)]
        no_reduction: bool,
    },

    /// Test whether a file can be fully and correctly compressed
    /// (for debugging only)
    #[command(arg_required_else_help = true)]
    Test {
        /// Input GFA file
        #[arg(required = true)]
        file: PathBuf,
    },

    /// Check whether a compressed file contains any duplicates
    /// (for debugging only)
    #[command(arg_required_else_help = true)]
    CheckCompressibility {
        /// Input GFA file
        #[arg(required = true)]
        file: PathBuf,
    },

    /// Reads in a file that already contains a grammar and tries to
    /// compress the file according to this grammar. This can be used
    /// to add new walks to an existing compressed graph.
    #[command(arg_required_else_help = true)]
    CompressAdditional {
        /// Input GFA file
        #[arg(required = true)]
        file: PathBuf,
    },
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Usage {
    Rule(UndirectedNodeId, usize),
    Haplotype(usize, usize),
    None,
    Multiple,
}

fn reduce_using_rule_utility(
    grammar: &mut DeterministicHashMap<UndirectedNodeId, Vec<NodeId>>,
    haplotypes: &mut [Vec<NodeId>],
) {
    let mut usages =
        DeterministicHashMap::with_capacity_and_hasher(grammar.len(), BuildHasherDefault::new());
    for k in grammar.keys() {
        usages.insert(*k, Usage::None);
    }

    for (key, rule) in grammar.iter() {
        for (idx, node) in rule.iter().enumerate() {
            if node.is_meta_node() {
                if usages[&node.0] == Usage::None {
                    usages
                        .entry(node.0)
                        .and_modify(|e| *e = Usage::Rule(*key, idx));
                } else {
                    usages.entry(node.0).and_modify(|e| *e = Usage::Multiple);
                }
            }
        }
    }

    for (hap_idx, haplo) in haplotypes.iter().enumerate() {
        for (idx, node) in haplo.iter().enumerate() {
            if node.is_meta_node() {
                if usages[&node.0] == Usage::None {
                    usages
                        .entry(node.0)
                        .and_modify(|e| *e = Usage::Haplotype(hap_idx, idx));
                } else {
                    usages.entry(node.0).and_modify(|e| *e = Usage::Multiple);
                }
            }
        }
    }

    // Get the necessary edits for haplotypes and walks
    let mut hap_edits: HashMap<usize, Vec<(usize, UndirectedNodeId)>> = HashMap::new();
    let mut rule_edits: HashMap<UndirectedNodeId, Vec<(usize, UndirectedNodeId)>> = HashMap::new();
    let mut rules_to_delete: Vec<UndirectedNodeId> = Vec::new();
    for (key, usage) in usages {
        match usage {
            Usage::Haplotype(hap_idx, idx) => {
                hap_edits.entry(hap_idx).or_default().push((idx, key));
                rules_to_delete.push(key);
            }
            Usage::Rule(rule_idx, idx) => {
                rule_edits.entry(rule_idx).or_default().push((idx, key));
                rules_to_delete.push(key);
            }
            _ => {}
        }
    }

    // Sort both edits descending to avoid messing up the indices when replacing
    for val in hap_edits.values_mut() {
        val.sort_by(|a, b| b.cmp(a));
    }
    for val in rule_edits.values_mut() {
        val.sort_by(|a, b| b.cmp(a));
    }

    // Do rule edits
    let mut queue = VecDeque::from(rule_edits.keys().copied().collect_vec());
    while let Some(curr) = queue.pop_front() {
        // If any of the elements in this rule will be replaced later, skip it for now
        if grammar[&curr].iter().any(|n| rule_edits.contains_key(&n.0)) {
            queue.push_back(curr);
        } else {
            // Directly remove the rule edit, we don't need this later
            let edit = rule_edits.remove(&curr).unwrap();
            for (i, _) in edit {
                let item = grammar[&curr][i];
                let rule_text = grammar[&item.0].clone();
                if item.is_forward() {
                    grammar.get_mut(&curr).unwrap().splice(i..=i, rule_text);
                } else {
                    grammar.get_mut(&curr).unwrap().splice(
                        i..=i,
                        rule_text.into_iter().rev().map(|n| n.flip()).collect_vec(),
                    );
                }
            }
        }
    }

    // Do haplotype edits
    for (key, edit) in hap_edits {
        for (i, _) in edit {
            let item = haplotypes[key][i];
            let rule_text = &grammar[&item.0];
            if item.is_forward() {
                haplotypes[key].splice(i..=i, rule_text.clone());
            } else {
                haplotypes[key].splice(
                    i..=i,
                    rule_text.iter().rev().map(|n| n.flip()).collect_vec(),
                );
            }
        }
    }

    // Delete removed rules
    for rule in rules_to_delete {
        grammar.remove(&rule);
    }
}

fn main() -> Result<()> {
    env_logger::Builder::from_env(Env::default().default_filter_or("info")).init();
    let args = Cli::parse();

    match args.command {
        Commands::CompressFull {
            file,
            prefix: _,
            no_reduction,
        } => {
            log::info!("Compressing file {:?}", file);
            let (haplotypes, mut node_registry, singleton_haplotypes) =
                parse_file_to_digrams(&file, true)?;
            let (haplotype_names, haplotypes): (Vec<PathSegment>, Vec<Vec<LocalizedDigram>>) =
                haplotypes.into_iter().unzip();
            let number_of_paths = haplotypes.len();
            log::info!("Parsed {} haplotypes", number_of_paths);
            let mut d = DigramOccurrences::from(haplotypes);
            let grammar = build_grammar(&mut d, &mut node_registry);
            d.log_sizes();
            let rev_reg: ReverseNodeRegistry = node_registry.into();
            let mut haplotype_walks =
                get_haplotype_walks(&d, &grammar, &singleton_haplotypes, number_of_paths);
            let mut grammar: DeterministicHashMap<UndirectedNodeId, Vec<NodeId>> = grammar
                .into_iter()
                .map(|(k, v0, v1, _)| (k.0, vec![v0, v1]))
                .collect();
            if !no_reduction {
                let grammar_len = grammar.len();
                reduce_using_rule_utility(&mut grammar, &mut haplotype_walks);
                log::info!(
                    "Reduced grammar by {} rules ({}%)",
                    grammar_len - grammar.len(),
                    (grammar_len - grammar.len()) as f64 / grammar_len as f64 * 100.0
                );
            }
            print_grammar_simple(&grammar, &rev_reg);
            let named_haplotypes: Vec<(PathSegment, Vec<NodeId>)> =
                haplotype_names.into_iter().zip(haplotype_walks).collect();
            print_walks(&named_haplotypes, &rev_reg);
        }
        Commands::Decompress { file, use_p_lines } => {
            log::info!("Decompressing file {:?}", file);
            let (haplotypes, node_registry, grammar) =
                parse_file_to_haplotypes_with_grammar(&file, true)?;
            let rev_reg: ReverseNodeRegistry = node_registry.into();
            decode_and_print_walks(haplotypes, &grammar, &rev_reg, use_p_lines);
        }
        Commands::CompressPartial {
            file,
            prefix: _,
            ratio,
            no_reduction,
        } => {
            log::info!("Compressing file {:?}", file);
            let (grammar, compressed_paths, node_registry, rev_reg) = {
                let (haplotypes, mut node_registry, singleton_haplotypes, compressed_paths) =
                    parse_file_to_digrams_ratio_based(&file, true, ratio)?;
                let (haplotype_names, haplotypes): (Vec<PathSegment>, Vec<Vec<LocalizedDigram>>) =
                    haplotypes.into_iter().unzip();
                let number_of_paths = haplotypes.len();
                log::info!("Parsed {} haplotypes", number_of_paths);
                let mut d = DigramOccurrences::from(haplotypes);
                let grammar = build_grammar(&mut d, &mut node_registry);
                let mut haplotype_walks =
                    get_haplotype_walks(&d, &grammar, &singleton_haplotypes, number_of_paths);
                let mut grammar: DeterministicHashMap<UndirectedNodeId, Vec<NodeId>> = grammar
                    .into_iter()
                    .map(|(k, v0, v1, _)| (k.0, vec![v0, v1]))
                    .collect();
                let rev_reg: ReverseNodeRegistry = node_registry.clone().into();
                if !no_reduction {
                    let grammar_len = grammar.len();
                    reduce_using_rule_utility(&mut grammar, &mut haplotype_walks);
                    log::info!(
                        "Reduced grammar by {} rules ({}%)",
                        grammar_len - grammar.len(),
                        (grammar_len - grammar.len()) as f64 / grammar_len as f64 * 100.0
                    );
                }
                print_grammar_simple(&grammar, &rev_reg);

                let named_haplotypes: Vec<(PathSegment, Vec<NodeId>)> =
                    haplotype_names.into_iter().zip(haplotype_walks).collect();
                print_walks(&named_haplotypes, &rev_reg);
                (grammar, compressed_paths, node_registry, rev_reg)
            };

            log::info!(
                "Created registries of size {} and {}",
                node_registry.deep_size_of() as f64 / 1e9,
                rev_reg.deep_size_of() as f64 / 1e9
            );
            let grammar = grammar
                .into_iter()
                .map(|(k, v)| (NodeId(k, helpers::utils::Orientation::Forward), v))
                .collect();
            compress_remaining_file(&file, grammar, &compressed_paths, &node_registry, &rev_reg)?;
        }
        Commands::CompressAdditional { file } => {
            log::info!("Reading file {:?}", file);
            let (haplotypes, node_registry, grammar) =
                parse_file_to_haplotypes_with_grammar(&file, true)?;
            let rev_reg: ReverseNodeRegistry = node_registry.into();
            let grammar: DeterministicHashMap<UndirectedNodeId, Vec<NodeId>> = grammar
                .into_iter()
                .map(|(k, (a, b))| (k, vec![a, b]))
                .collect();
            print_grammar_simple(&grammar, &rev_reg);
            let grammar = grammar
                .into_iter()
                .map(|(k, v)| (NodeId(k, helpers::utils::Orientation::Forward), v))
                .collect();
            let named_walks = compress(haplotypes, grammar);
            print_walks(&named_walks, &rev_reg);
        }
        Commands::Test { file } => {
            log::info!("Testing on file {:?}", file);
            let (haplotypes, mut node_registry, singleton_haplotypes) =
                parse_file_to_digrams(&file, false)?;
            let (_, haplotypes): (Vec<PathSegment>, Vec<Vec<LocalizedDigram>>) =
                haplotypes.into_iter().unzip();
            let number_of_paths = haplotypes.len();
            log::info!("Parsed {} haplotypes", number_of_paths);
            let mut d = DigramOccurrences::from(haplotypes);
            log::info!("Constructed occurrences of {} digrams", d.total_len());
            let grammar = build_grammar(&mut d, &mut node_registry);
            log::info!("Built grammar with {} rules", grammar.len());

            let haplotype_walks =
                get_haplotype_walks(&d, &grammar, &singleton_haplotypes, number_of_paths);
            log::info!("Encoded {} haplotypes", haplotype_walks.len());

            let grammar: Grammar = grammar
                .into_iter()
                .map(|(name, u, v, _)| (name.get_undirected(), (u, v)))
                .collect();
            check_incompressibility(&haplotype_walks, &grammar, &node_registry);

            compare_file(&file, &haplotype_walks, &grammar, &node_registry);
        }
        Commands::CheckCompressibility { file } => {
            log::info!("Decompressing file {:?}", file);
            let (haplotypes, node_registry, grammar) =
                parse_file_to_haplotypes_with_grammar(&file, false)?;
            let unnamed_haplotypes = haplotypes.into_iter().map(|(_, w)| w).collect_vec();
            check_incompressibility(&unnamed_haplotypes, &grammar, &node_registry);
        }
    }
    Ok(())
}

fn check_incompressibility(walks: &[Vec<NodeId>], grammar: &Grammar, node_registry: &NodeRegistry) {
    let mut digrams: HashSet<CanonicalDigram> = HashSet::new();
    let mut total_seen_digrams = 0;
    let rev_node: ReverseNodeRegistry = node_registry.clone().into();

    for rule in grammar.values() {
        let digram: CanonicalDigram = Digram::new(rule.0, rule.1).into();
        if !digrams.insert(digram.clone()) {
            log::error!(
                "Digram {} (internal: {}) - {} (internal: {}) seen twice in grammar",
                rev_node.get_directed_name(digram.0),
                digram.0 .0 .0,
                rev_node.get_directed_name(digram.1),
                digram.1 .0 .0,
            );
        }
        total_seen_digrams += 1;
    }

    for walk in walks {
        for (u, v) in walk.iter().tuple_windows() {
            let digram: CanonicalDigram = Digram::new(*u, *v).into();
            if !digrams.insert(digram.clone()) {
                log::error!(
                    "Digram {} (internal: {}) - {} (internal: {}) seen twice",
                    rev_node.get_directed_name(digram.0),
                    digram.0 .0 .0,
                    rev_node.get_directed_name(digram.1),
                    digram.1 .0 .0,
                );
            }
            total_seen_digrams += 1;
        }
    }

    if digrams.len() != total_seen_digrams {
        log::error!(
            "{} digrams seen multiple times",
            total_seen_digrams - digrams.len()
        );
    } else {
        println!("Grammar and haplotypes cannot be compressed any further");
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn n(x: u32) -> NodeId {
        NodeId(
            UndirectedNodeId(x, false),
            helpers::utils::Orientation::Forward,
        )
    }

    fn m(x: u32) -> NodeId {
        NodeId(
            UndirectedNodeId(x, true),
            helpers::utils::Orientation::Forward,
        )
    }

    fn rn(x: u32) -> NodeId {
        NodeId(
            UndirectedNodeId(x, false),
            helpers::utils::Orientation::Backward,
        )
    }

    fn rm(x: u32) -> NodeId {
        NodeId(
            UndirectedNodeId(x, true),
            helpers::utils::Orientation::Backward,
        )
    }

    #[test]
    fn test_reduce_every_rule() {
        let mut grammar = DeterministicHashMap::default();
        grammar.insert(UndirectedNodeId(7, true), vec![n(1), n(2)]);
        grammar.insert(UndirectedNodeId(8, true), vec![m(7), n(3)]);
        grammar.insert(UndirectedNodeId(9, true), vec![m(8), n(4)]);
        let mut haplos = vec![vec![m(9), n(5), n(6)]];
        reduce_using_rule_utility(&mut grammar, &mut haplos);
        assert_eq!(grammar.len(), 0);
        assert_eq!(haplos, vec![vec![n(1), n(2), n(3), n(4), n(5), n(6)]]);
    }

    #[test]
    fn test_reduce_one_rule() {
        let mut grammar = DeterministicHashMap::default();
        grammar.insert(UndirectedNodeId(6, true), vec![n(1), n(2)]);
        grammar.insert(UndirectedNodeId(7, true), vec![m(6), n(3)]);
        let mut haplos = vec![vec![m(7)], vec![m(7)]];
        reduce_using_rule_utility(&mut grammar, &mut haplos);
        assert_eq!(haplos, vec![vec![m(7)], vec![m(7)]]);
        assert_eq!(grammar.len(), 1);
        assert_eq!(grammar[&UndirectedNodeId(7, true)], vec![n(1), n(2), n(3)]);
    }

    #[test]
    fn test_reduce_one_rule_rev() {
        let mut grammar = DeterministicHashMap::default();
        grammar.insert(UndirectedNodeId(6, true), vec![n(1), n(2)]);
        grammar.insert(UndirectedNodeId(7, true), vec![rm(6), n(3)]);
        let mut haplos = vec![vec![m(7)], vec![rm(7)]];
        reduce_using_rule_utility(&mut grammar, &mut haplos);
        assert_eq!(haplos, vec![vec![m(7)], vec![rm(7)]]);
        assert_eq!(grammar.len(), 1);
        assert_eq!(
            grammar[&UndirectedNodeId(7, true)],
            vec![rn(2), rn(1), n(3)]
        );
    }
}
