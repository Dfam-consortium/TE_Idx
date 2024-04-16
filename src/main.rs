use std::num::ParseIntError;
use std::path::Path;

use clap::builder::PossibleValuesParser;
use clap::{Parser, Subcommand};

use te_idx::bgzf_filter;
use te_idx::idx_query;
use te_idx::json_query;
use te_idx::prep_beds;
use te_idx::process_json;
use te_idx::read_family_assembly_annotations;
mod idx;

use te_idx::{ASSEMBLY_DIR, BENCHMARK_DIR, DATA_DIR, MASKS_DIR, MOD_LEN_DIR, SEQUENCE_DIR};

const INDEX_DATA_TYPES: [&str; 3] = [ASSEMBLY_DIR, BENCHMARK_DIR, MASKS_DIR];
const JSON_DATA_TYPES: [&str; 2] = [MOD_LEN_DIR, SEQUENCE_DIR];

#[derive(Parser)]
#[command(author, version, about)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Option<Commands>,
}

#[derive(Subcommand)]
pub enum Commands {
    /// Search an individual .bed.bgz file for matches in a specified column
    BgzfFilter {
        /// Path to .bed.bgz file to be searched. Assumed to be CSV/TSV
        #[arg(long, short, verbatim_doc_comment)]
        infile: String,
        /// Column number to be searched. 1-indexed
        #[arg(long, short, verbatim_doc_comment)]
        position: usize,
        /// Term to be searched for. If absent, all rows will be returned
        #[arg(long, short, verbatim_doc_comment)]
        term: Option<String>,
        /// Path to file to save filtered data. Should end in .bed.bgz
        #[arg(long, short, verbatim_doc_comment)]
        outfile: Option<String>,
    },
    /// Build and/or Search Index file for grouped .bed.bgz files
    Idx {
        /// Name of assembly/assembly folder
        #[arg(short, long, verbatim_doc_comment)]
        assembly: String,
        /// Type of data to be indexed/searched
        #[arg(short, long, verbatim_doc_comment)]
        #[clap(value_parser = PossibleValuesParser::new(INDEX_DATA_TYPES), required(true))]
        data_type: String,
        /// Trigger index build
        #[arg(short, long, verbatim_doc_comment)]
        build: bool,
        /// Search index. Must take the form: Chrom-Start-End
        #[arg(short, long, verbatim_doc_comment)]
        search: Option<String>,
        /// Optional: Only return hits matching accession
        #[arg(short, long, verbatim_doc_comment)]
        family: Option<String>,
        /// Optional: Only return NRPH hits
        #[arg(short, long, verbatim_doc_comment)]
        nrph: bool,
    },
    /// Split TSV files into compressed BED files by accession
    PrepBeds {
        /// Name of assembly/assembly folder
        #[arg(short, long, verbatim_doc_comment)]
        assembly: String,
        /// Input file from buildFullRegion.py. Should be in accession order
        #[arg(short, long, verbatim_doc_comment)]
        in_tsv: String,
        /// Type of data to be prepped
        #[arg(short, long, verbatim_doc_comment)]
        #[clap(value_parser = PossibleValuesParser::new(INDEX_DATA_TYPES), required(true))]
        data_type: String,
    },
    /// Convert JSON exports from PHPMyAdmin to a useable JSON file
    ProcessJSON {
        /// Input JSON file, exported from PHPMyAdmin
        #[arg(short, long, verbatim_doc_comment)]
        in_file: String,
        /// Key attribute, such as accession, to be used in the new map
        #[arg(short, long, verbatim_doc_comment)]
        key: String,
        /// Optional: Output file. If not provided, will print to stdout
        #[arg(long, short)]
        outfile: Option<String>,
    },
    /// Search indexed BED files for all hits within a range
    IdxQuery {
        /// Name of assembly/assembly folder
        #[arg(short, long, verbatim_doc_comment)]
        assembly: String,
        /// Type of data to be searched
        #[arg(short, long, verbatim_doc_comment)]
        #[clap(value_parser = PossibleValuesParser::new(INDEX_DATA_TYPES), required(true))]
        data_type: String,
        /// chromosome number/accession
        #[arg(short, long, verbatim_doc_comment)]
        chrom: String,
        /// start position
        #[arg(short, long, verbatim_doc_comment)]
        start: u64,
        /// end position
        #[arg(short, long, verbatim_doc_comment)]
        end: u64,
        /// Optional: Only return hits matching accession
        #[arg(short, long, verbatim_doc_comment)]
        family: Option<String>,
        /// Only return NRPH hits
        #[arg(short, long, verbatim_doc_comment)]
        nrph: bool,
    },
    /// Retrieve information from a processed JSON file
    JsonQuery {
        /// Name of assembly/assembly folder
        #[arg(short, long, verbatim_doc_comment)]
        assembly: String,
        /// Type of data to be searched
        #[arg(short, long, verbatim_doc_comment)]
        #[clap(value_parser = PossibleValuesParser::new(JSON_DATA_TYPES), required(true))]
        data_type: String,
        /// Key value to search by
        #[arg(short, long, verbatim_doc_comment)]
        key: String,
        /// Target Attribute to return
        #[arg(short, long, verbatim_doc_comment)]
        target: String,
    },
    /// Read all or NRPH only family annotations for an assembly
    ReadFamilyAssemblyAnnotations {
        /// Family Accession
        #[arg(short, long, verbatim_doc_comment)]
        id: String,
        /// Name of assembly/assembly folder
        #[arg(short, long, verbatim_doc_comment)]
        assembly_id: String,
        /// Only Return NRPH hits
        #[arg(short, long, verbatim_doc_comment)]
        nrph: bool,
        /// Optional: Output file
        #[arg(long, short)]
        outfile: Option<String>,
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
    let cli = Cli::parse();

    match &cli.command {
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
        Some(Commands::IdxQuery {
            assembly,
            data_type,
            chrom,
            start,
            end,
            family,
            nrph,
        }) => idx_query(assembly, data_type, chrom, *start, *end, family, nrph),
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
        Some(Commands::ProcessJSON {
            in_file,
            key,
            outfile,
        }) => process_json(in_file, key, outfile).expect("JSON Parse Failed"),
        None => {}
    }
}
