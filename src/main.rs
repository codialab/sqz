use anyhow::Result;
use clap::{Parser, Subcommand};
use env_logger::Env;
use helpers::{digram_occurrences::DigramOccurrences, utils::LocalizedDigram};
use parser::parse_file_to_haplotypes;
use std::{collections::HashMap, path::PathBuf};

use crate::{
    encoding::get_haplotype_walks,
    grammar_building::{Rule, build_grammar},
    helpers::{PathSegment, ReverseNodeRegistry, utils::NodeId},
    printing::{print_grammar, print_walks},
};

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

    #[command(arg_required_else_help = true)]
    Test {
        #[arg(required = true)]
        file: String,
    },
}

pub fn get_digrams(haplotypes: &[Vec<NodeId>]) -> (Vec<Vec<LocalizedDigram>>, HashMap<usize, NodeId>) {
    let mut singleton_haplotypes = HashMap::new();
    let mut digrams = vec![Vec::new(); haplotypes.len()];

    for (path_idx, haplotype) in haplotypes.iter().enumerate() {
        if haplotype.len() == 1 {
            singleton_haplotypes.insert(path_idx, haplotype[0]);
            continue;
        }
        let mut address_counter = 0;
        for 
    }

    (digrams, singleton_haplotypes)
}

fn main() -> Result<()> {
    env_logger::Builder::from_env(Env::default().default_filter_or("info")).init();
    let args = Cli::parse();

    match args.command {
        Commands::Compress { file, prefix: _ } => {
            log::info!("Compressing file {:?}", file);
            let (haplotypes, mut node_registry) =
                parse_file_to_haplotypes(&file, true)?;
            let (haplotype_names, haplotypes): (Vec<PathSegment>, Vec<Vec<NodeId>>) =
                haplotypes.into_iter().unzip();
            let (digrams, singleton_haplotypes) = get_digrams(haplotypes);
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
            use_p_lines: _,
        } => {
            log::info!("Decompressing file {}", file);
        }
        Commands::Test { file } => {
            log::info!("Testing on file {}", file);
        }
    }
    Ok(())
}
