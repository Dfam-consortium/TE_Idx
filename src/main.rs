use std::num::ParseIntError;
use std::path::Path;

use clap::builder::PossibleValuesParser;
use clap::{Parser, Subcommand};

use te_idx::bgzf_filter;
use te_idx::json_query;
use te_idx::nhmmer_query;
use te_idx::prep_beds;
use te_idx::process_json;
use te_idx::read_family_assembly_annotations;
use te_idx::trf_query;

mod idx;

use te_idx::{ASSEMBLY_DIR, BENCHMARK_DIR, DATA_DIR, MASKS_DIR};

const INDEX_DATA_TYPES: [&str; 3] = [ASSEMBLY_DIR, BENCHMARK_DIR, MASKS_DIR];

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
        term: Option<String>,
        #[arg(long, short)]
        outfile: Option<String>,
    },
    Idx {
        #[arg(short, long)]
        assembly: String,

        #[arg(short, long)]
        #[clap(value_parser = PossibleValuesParser::new(INDEX_DATA_TYPES), required(true))]
        data_type: String,

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
    PrepBeds {
        #[arg(short, long)]
        assembly: String,

        #[arg(short, long)]
        in_tsv: String,

        #[arg(short, long)]
        #[clap(value_parser = PossibleValuesParser::new(INDEX_DATA_TYPES), required(true))]
        data_type: String,
    },
    NhmmerQuery {
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
    TrfQuery {
        #[arg(short, long)]
        assembly: String,
        #[arg(short, long)]
        chrom: String,
        #[arg(short, long)]
        start: u64,
        #[arg(short, long)]
        end: u64,
    },
    JsonQuery {
        #[arg(short, long)]
        assembly: String,

        #[arg(short, long)]
        data_type: String,

        #[arg(short, long)]
        key: String,

        #[arg(short, long)]
        target: String,
    },
    ReadFamilyAssemblyAnnotations {
        #[arg(short, long)]
        id: String,
        #[arg(short, long)]
        assembly_id: String,
        #[arg(short, long)]
        nrph: bool,
        #[arg(long, short)]
        outfile: Option<String>,
    },
    ProcessJSON {
        #[arg(short, long)]
        in_file: String,

        #[arg(short, long)]
        key: String,
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
            data_type,
            build,
            search,
            family,
            nrph,
        }) => {
            let (filenames, bgz_dir, mut contig_index, index_file) =
                match idx::prep_idx(&format!("{}/{}", &DATA_DIR, &assembly), data_type) {
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
                let results = idx::search_idx(
                    &filenames,
                    &bgz_dir,
                    &mut contig_index,
                    &index_file,
                    &q_contig.to_string(),
                    start,
                    end,
                    family,
                    *nrph,
                    false,
                );
                println!("{:?}", results);
            }
        }
        Some(Commands::PrepBeds {
            assembly,
            in_tsv,
            data_type,
        }) => match prep_beds(assembly, in_tsv, data_type) {
            Ok(()) => println!("Bed Files Created - {}", data_type),
            Err(e) => panic!("{:?}", e),
        },
        Some(Commands::NhmmerQuery {
            assembly,
            chrom,
            start,
            end,
            family,
            nrph,
        }) => nhmmer_query(assembly, chrom, *start, *end, family, nrph),
        Some(Commands::TrfQuery {
            assembly,
            chrom,
            start,
            end,
        }) => trf_query(assembly, chrom, *start, *end).expect("Mask Read Failed"),
        Some(Commands::JsonQuery {
            assembly,
            data_type,
            key,
            target,
        }) => {
            let ans = json_query(assembly, data_type, key, target).expect("JSON Read Failed");
            println!("{}", ans)
        }
        Some(Commands::ReadFamilyAssemblyAnnotations {
            id,
            assembly_id,
            nrph,
            outfile,
        }) => read_family_assembly_annotations(id, assembly_id, nrph, outfile),
        Some(Commands::ProcessJSON { in_file, key }) => {
            process_json(in_file, key).expect("JSON Parse Failed")
        }
        None => {}
    }
}
