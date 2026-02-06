use anyhow::Result;
use clap::{Parser, Subcommand};
use deepsize::DeepSizeOf;
use env_logger::Env;
use helpers::{digram_occurrences::DigramOccurrences, utils::LocalizedDigram};
use itertools::Itertools;
use std::{collections::HashSet, path::PathBuf};

use crate::{
    compressing::compress_remaining_file,
    decoding::decode_and_print_walks,
    encoding::get_haplotype_walks,
    grammar_building::{build_grammar, Rule},
    helpers::{
        utils::{CanonicalDigram, Digram, NodeId},
        DeterministicHashMap, NodeRegistry, PathSegment, ReverseNodeRegistry,
    },
    parser::{
        compare_file, parse_file_to_digrams, parse_file_to_digrams_ratio_based,
        parse_file_to_haplotypes_with_grammar, Grammar,
    },
    printing::{print_grammar, print_grammar_simple, print_walks},
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
    },

    #[command(arg_required_else_help = true)]
    Test {
        #[arg(required = true)]
        file: PathBuf,
    },
}

fn main() -> Result<()> {
    env_logger::Builder::from_env(Env::default().default_filter_or("info")).init();
    let args = Cli::parse();

    match args.command {
        Commands::CompressFull { file, prefix: _ } => {
            log::info!("Compressing file {:?}", file);
            let (haplotypes, mut node_registry, singleton_haplotypes) =
                parse_file_to_digrams(&file, true)?;
            let (haplotype_names, haplotypes): (Vec<PathSegment>, Vec<Vec<LocalizedDigram>>) =
                haplotypes.into_iter().unzip();
            let number_of_paths = haplotypes.len();
            log::info!("Parsed {} haplotypes", number_of_paths);
            let mut d = DigramOccurrences::from(haplotypes);
            let grammar = build_grammar(&mut d, &mut node_registry);
            let rev_reg: ReverseNodeRegistry = node_registry.into();
            for (digram, occurrences) in &d {
                log::info!(
                    "| {}{} -> {:?}",
                    rev_reg.get_directed_name(digram.0),
                    rev_reg.get_directed_name(digram.1),
                    occurrences
                );
            }
            print_grammar(&grammar, &rev_reg, true);
            let haplotype_walks =
                get_haplotype_walks(&d, &grammar, &singleton_haplotypes, number_of_paths);
            print_walks(&haplotype_walks, &rev_reg, &haplotype_names);
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
                let haplotype_walks =
                    get_haplotype_walks(&d, &grammar, &singleton_haplotypes, number_of_paths);
                let grammar: DeterministicHashMap<NodeId, Vec<NodeId>> = grammar
                    .into_iter()
                    .map(|(k, v0, v1, _)| (k, vec![v0, v1]))
                    .collect();
                print_grammar_simple(&grammar);

                let rev_reg: ReverseNodeRegistry = node_registry.clone().into();
                print_walks(&haplotype_walks, &rev_reg, &haplotype_names);
                (grammar, compressed_paths, node_registry, rev_reg)
            };

            log::info!(
                "Created registries of size {} and {}",
                node_registry.deep_size_of() as f64 / 1e9,
                rev_reg.deep_size_of() as f64 / 1e9
            );
            compress_remaining_file(&file, grammar, &compressed_paths, &node_registry, &rev_reg)?;
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
                "Digram {:?} - {:?} seen twice in grammar",
                rev_node.get_directed_name(digram.0),
                rev_node.get_directed_name(digram.1)
            );
        }
        total_seen_digrams += 1;
    }

    for walk in walks {
        for (u, v) in walk.iter().tuple_windows() {
            let digram: CanonicalDigram = Digram::new(*u, *v).into();
            if !digrams.insert(digram.clone()) {
                log::error!(
                    "Digram {:?} - {:?} seen twice",
                    rev_node.get_directed_name(digram.0),
                    rev_node.get_directed_name(digram.1)
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
