use anyhow::Result;
use clap::{Parser, Subcommand};
use env_logger::Env;
use helpers::{digram_occurrences::DigramOccurrences, utils::LocalizedDigram};
use itertools::Itertools;
use std::{collections::HashSet, path::PathBuf};

use crate::{
    decoding::decode_and_print_walks, encoding::get_haplotype_walks, grammar_building::{Rule, build_grammar}, helpers::{PathSegment, ReverseNodeRegistry, utils::{CanonicalDigram, Digram, NodeId}}, parser::{Grammar, compare_file, parse_file_to_digrams, parse_file_to_haplotypes_with_grammar}, printing::{print_grammar, print_walks}
};

mod encoding;
mod decoding;
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
        file: PathBuf,

        #[arg(short = 'p', long)]
        use_p_lines: bool,
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
        Commands::Compress { file, prefix: _ } => {
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
        Commands::Decompress {
            file,
            use_p_lines,
        } => {
            log::info!("Decompressing file {:?}", file);
            let (haplotypes, node_registry, grammar) = parse_file_to_haplotypes_with_grammar(&file, true)?;
            let rev_reg: ReverseNodeRegistry = node_registry.into();
            decode_and_print_walks(haplotypes, &grammar, &rev_reg, use_p_lines);
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
            let grammar = build_grammar(&mut d, &mut node_registry);
            log::info!("Built grammar with {} rules", grammar.len());

            let haplotype_walks =
                get_haplotype_walks(&d, &grammar, &singleton_haplotypes, number_of_paths);
            log::info!("Encoded {} haplotypes", haplotype_walks.len());

            let grammar: Grammar = grammar.into_iter().map(|(name, u, v, _)| (name.get_undirected(), (u, v))).collect(); 
            check_incompressibility(&haplotype_walks, &grammar);

            compare_file(&file, &haplotype_walks, &grammar, &node_registry);
        }
    }
    Ok(())
}

fn check_incompressibility(walks: &[Vec<NodeId>], grammar: &Grammar) {
    let mut digrams: HashSet<CanonicalDigram> = HashSet::new();
    let mut total_seen_digrams = 0;

    for (_, rule) in grammar {
        let digram: CanonicalDigram = Digram::new(rule.0, rule.1).into();
        digrams.insert(digram);
        total_seen_digrams += 1;
    }

    for walk in walks {
        for (u, v) in walk.iter().tuple_windows() {
            let digram: CanonicalDigram = Digram::new(*u, *v).into();
            digrams.insert(digram);
            total_seen_digrams += 1;
        }
    }

    assert_eq!(digrams.len(), total_seen_digrams);
    println!("Grammar and haplotypes cannot be compressed any further");
}
