use anyhow::Result;
use clap::{Parser, Subcommand};
use helpers::{Freq, Haplotype, NeighborList, D};
use parser::parse_file_to_haplotypes;
use std::path::PathBuf;

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

fn get_data_structures(haplotypes: Vec<Haplotype>) -> (D, Freq, NeighborList, NeighborList) {
    let mut d = D::new();
    let mut neighbor_left = NeighborList::new();
    let mut neighbor_right = NeighborList::new();

    for (i, haplotype) in haplotypes.into_iter().enumerate() {
        for local_digram in haplotype {
            let (digram, address) = local_digram.split_to_canonical();
            d.insert(&digram, (i, address.clone()));
            neighbor_right.insert(
                (digram.get_u(), i, address.get_first()),
                (digram.get_v(), address.get_second()),
            );
            neighbor_left.insert(
                (digram.get_v(), i, address.get_second()),
                (digram.get_u(), address.get_first()),
            );
        }
    }

    let freq = d.get_freq();

    (d, freq, neighbor_left, neighbor_right)
}

fn main() -> Result<()> {
    env_logger::init();
    let args = Cli::parse();

    match args.command {
        Commands::Compress { file, prefix } => {
            log::info!("Compressing file {:?}", file);
            let haplotypes = parse_file_to_haplotypes(&file)?;
            log::info!("Parsed {} haplotypes", haplotypes.len());
            let (d, freq, neighbor_left, neighbor_right) = get_data_structures(haplotypes);
        }
        Commands::Decompress { file, use_p_lines } => {
            log::info!("Decompressing file {}", file);
        }
    }
    Ok(())
}
