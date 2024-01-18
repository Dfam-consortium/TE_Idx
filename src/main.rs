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
        #[arg(short, long)]
        proj_dir: String,

        // Build the index from the current hard-coded filename
        #[arg(short, long)]
        build: bool,

        // Search the index using a hard coded search range
        #[arg(short, long)]
        search: bool,
    },
}

fn main() {
    let cli = Cli::parse(); // TODO update template
 // Using builder interface to support a custom help template
    //     let cli = Command::new("te_idx").help_template(
    //         "
    // {before-help}{name} {version}
    // {author-with-newline}{about-with-newline}
    // Te_idx is a tool for indexing/searching Transposable Element annotation datasets.
    // ...
    // The index is based on the IGD linear binning approach with extensions for accessing
    // metadata from compressed BED files (blocked gzip).

    // {usage-heading} {usage}

    //   Examples:

    //     # Build the index
    //     ./te_idx --build

    //     # Search the index
    //     ./te_idx --search

    // {all-args}{after-help}
    // ",
    //     );

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
        Some(Commands::Idx {
            proj_dir,
            build,
            search,
        }) => {
            let (filenames, bgz_dir, contig_index, index_file) = idx::prep_idx(proj_dir);
            if *build {
                idx::build_idx(filenames, bgz_dir, contig_index, index_file)
            } else if *search {
                idx::search_idx(filenames, bgz_dir, contig_index, index_file)
            }
        }
        None => {}
    }
}
