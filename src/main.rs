use clap::{Parser, Subcommand};

use te_idx::bgzf_filter;
use te_idx::find_family;

mod idx;

#[derive(Parser)]
#[command(author, version, about, long_about = None)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Option<Commands>,
}

#[derive(Subcommand)]
pub enum Commands {
    BgzfFilter {
        #[arg(long, short)]
        infile: String,
        #[arg(long, short)]
        position: usize,
        #[arg(long, short)]
        term: String,
        #[arg(long, short)]
        outfile: String,
    },
    FindFamily {
        #[arg(long, short)]
        id: String,
        #[arg(long, short)]
        assembly: String,
    },
    Idx {
        // Build the index from the current hard-coded filename
        #[arg(short, long)]
        build: bool,

        // Search the index using a hard coded search range
        #[arg(short, long)]
        search: bool,
    }
}

fn main() {
    let cli = Cli::parse();

    match &cli.command {
        Some(Commands::FindFamily { id, assembly }) => find_family(id, assembly),
        Some(Commands::BgzfFilter {
            infile,
            position,
            term,
            outfile,
        }) => {
            bgzf_filter(infile, position, term, outfile).expect("Filter Failed");
        }
        None => {}
    }
}
