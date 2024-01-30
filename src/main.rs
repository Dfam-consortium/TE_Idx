use std::num::ParseIntError;
use std::path::Path;

use clap::{Parser, Subcommand};

use te_idx::bgzf_filter;
// use te_idx::find_family;
use te_idx::read_annotations;
use te_idx::read_family_assembly_annotations;

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
        // term: String,
        term: Option<String>,
        #[arg(long, short)]
        // outfile: String,
        outfile: Option<String>,
    },
    // FindFamily {
    //     #[arg(long, short)]
    //     id: String,
    //     #[arg(long, short)]
    //     assembly: String,
    // },
    Idx {
        #[arg(short, long)]
        assembly: String,

        // Boolean flag to build/rebuild index
        #[arg(short, long)]
        build: bool,

        // Search the index using a search range
        #[arg(short, long)]
        search: Option<String>,

        #[arg(short, long)]
        family: Option<String>,

        #[arg(short, long)]
        nrph: bool,
    },
    ReadAnnotations {
        #[arg(short, long)]
        assembly: String,
        #[arg(short, long)]
        chrom: String,
        #[arg(short, long)]
        start: u64,
        #[arg(short, long)]
        end: u64,
        #[arg(short, long)]
        family: Option<String>,
        #[arg(short, long)]
        nrph: bool,
    },
    ReadFamilyAssemblyAnnotations {
        #[arg(short, long)]
        id: String,
        #[arg(short, long)]
        assembly_id: String,
        #[arg(short, long)]
        nrph: bool,
    },
}

fn parse_coordinates(coord_str: &str) -> Result<(u64, u64, u64), ParseIntError> {
    let coords: Vec<u64> = coord_str
        .split('-')
        .map(|e| match e.parse::<u64>() {
            Ok(val) => val,
            Err(e) => panic!("ERROR: Invalid Coordinates - {:?}", e),
        })
        .collect();
    if coords.len() != 3 {
        panic!("ERROR: Invalid Coordinates: Need to be in the form of chrom-start-end")
    } else {
        return Ok((coords[0], coords[1], coords[2]));
    }
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
        // Some(Commands::FindFamily { id, assembly }) => find_family(id, assembly),
        Some(Commands::BgzfFilter {
            infile,
            position,
            term,
            outfile,
        }) => {
            bgzf_filter(infile, position, term, outfile).expect("Filter Failed");
        }
        Some(Commands::Idx {
            assembly,
            build,
            search,
            family,
            nrph,
        }) => {
            let (filenames, bgz_dir, mut contig_index, index_file) = match idx::prep_idx(assembly) {
                Ok(res) => res,
                Err(e) => panic!("Search Prep Failed, Index may not exist - {:?}", e),
            };
            if *build {
                idx::build_idx(&filenames, &bgz_dir, &mut contig_index, &index_file)
            }
            if search.is_some() {
                if !Path::new(&index_file).exists() {
                    panic!(
                        "Index {:?} does not exist. Run the build command first",
                        &index_file
                    )
                }
                let (q_contig, start, end) = match parse_coordinates(search.as_deref().unwrap()) {
                    Ok(res) => res,
                    Err(e) => panic!("{:?}", e),
                };
                let _results = idx::search_idx(
                    &filenames,
                    &bgz_dir,
                    &mut contig_index,
                    &index_file,
                    &q_contig.to_string(),
                    start,
                    end,
                    family,
                    *nrph,
                );
            }
        }
        Some(Commands::ReadAnnotations {
            assembly,
            chrom,
            start,
            end,
            family,
            nrph,
        }) => read_annotations(assembly, chrom, *start, *end, family, nrph),
        Some(Commands::ReadFamilyAssemblyAnnotations {
            id,
            assembly_id,
            nrph,
        }) => read_family_assembly_annotations(id, assembly_id, nrph),
        None => {}
    }
}
