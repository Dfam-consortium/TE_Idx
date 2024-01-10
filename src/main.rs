use clap::{Parser, Subcommand};

use te_idx::bgzf_filter;
use te_idx::find_family;

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
        outfile: String
    },
    FindFamily {
        #[arg(long, short)]
        id: String,
    }
}

fn main() {
    let cli = Cli::parse();

    match &cli.command {
        Some(Commands::FindFamily {id}) => {
            find_family(id)
        },
        Some(Commands::BgzfFilter {infile, position, term, outfile}) => {
            bgzf_filter(infile, position, term, outfile).expect("Filter Failed");
        }
        None => {}
    }
}
