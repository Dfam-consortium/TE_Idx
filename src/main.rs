use clap::builder::PossibleValuesParser;
use clap::{Parser, Subcommand};

use te_idx::bgzf_filter;
use te_idx::idx_query;
use te_idx::json_query;
use te_idx::prep_beds;
use te_idx::prepare_assembly;
use te_idx::read_family_assembly_annotations;
use te_idx::get_chrom_id;
mod idx;

use te_idx::{DATA_DIR, INDEX_DATA_TYPES, JSON_DATA_TYPES};

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
        /// Name of assembly/assembly folder
        #[arg(short, long, verbatim_doc_comment)]
        assembly: String,
        /// Type of data to be indexed/searched
        #[arg(short, long, verbatim_doc_comment)]
        #[clap(value_parser = PossibleValuesParser::new(INDEX_DATA_TYPES), required(true))]
        data_type: String,
        /// Family name, corresponds to compressed TSV file prefix
        #[arg(long, short, verbatim_doc_comment)]
        fam: String,
        /// Column number to be searched. 1-indexed
        #[arg(long, short, verbatim_doc_comment)]
        position: usize,
        /// Term to be searched for. If absent, all rows will be returned
        #[arg(long, short, verbatim_doc_comment)]
        term: Option<String>,
        /// Path to file to save filtered data. Should end in .bed.bgz
        #[arg(long, short, verbatim_doc_comment)]
        outfile: Option<String>,
        /// Flag to reformat the feilds to match Dfam.org download file format
        #[arg(long, short, verbatim_doc_comment)]
        web_fmt: bool,
    },
    /// Build file for grouped .bed.bgz files
    BuildIdx {
        /// Name of assembly/assembly folder
        #[arg(short, long, verbatim_doc_comment)]
        assembly: String,
        /// Type of data to be indexed/searched
        #[arg(short, long, verbatim_doc_comment)]
        #[clap(value_parser = PossibleValuesParser::new(INDEX_DATA_TYPES), required(true))]
        data_type: String,
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
    /// Given an assembly name, check for and process all present exports
    PrepareAssembly {
        /// Name of assembly/assembly folder
        #[arg(short, long, verbatim_doc_comment)]
        assembly: String,
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
    GetChromID {
        /// Name of assembly/assembly folder
        #[arg(short, long, verbatim_doc_comment)]
        assembly: String,
        /// Query to get ID for
        #[arg(long, short)]
        query: String,
    },
}

fn main() {
    let cli = Cli::parse();

    match &cli.command {
        Some(Commands::BgzfFilter {
            assembly,
            data_type,
            fam,
            position,
            term,
            outfile,
            web_fmt,
        }) => {
            bgzf_filter(
                assembly, data_type, fam, position, term, outfile, *web_fmt, None,
            )
            .expect("Filter Failed");
        }
        Some(Commands::BuildIdx {
            assembly,
            data_type,
        }) => {
            let (filenames, bgz_dir, mut contig_index, index_file) =
                match idx::prep_idx(&format!("{}/{}", &DATA_DIR, &assembly), data_type) {
                    Ok(res) => res,
                    Err(e) => panic!(
                        "Search Prep Failed, Assembly or Data Type May Not Exist - {:?}",
                        e
                    ),
                };
            idx::build_idx(&filenames, &bgz_dir, &mut contig_index, &index_file)
                .expect("Indexing Failed")
        }
        Some(Commands::PrepBeds {
            assembly,
            in_tsv,
            data_type,
        }) => match prep_beds(assembly, in_tsv, data_type, None) {
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
        }) => {
            let result = idx_query(assembly, data_type, chrom, *start, *end, family, nrph, None)
                .expect("Index Query Failed");
            println!("{}", result)
        }
        Some(Commands::JsonQuery {
            assembly,
            data_type,
            key,
            target,
        }) => {
            let ans = json_query(assembly, data_type, key, target, None).expect("JSON Read Failed");
            println!("{}", ans)
        }
        Some(Commands::ReadFamilyAssemblyAnnotations {
            id,
            assembly_id,
            nrph,
            outfile,
        }) => {
            let _res = read_family_assembly_annotations(id, assembly_id, nrph, outfile, None);
        }
        Some(Commands::PrepareAssembly { assembly }) => prepare_assembly(assembly, None, None)
            .expect(format!("Assembly Prep for {} Failed", &assembly).as_str()),
        Some(Commands::GetChromID { assembly , query}) => {println!("{}", get_chrom_id(assembly, query));},
        None => {}
    }
}

// OLD CODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// use te_idx::process_json;
//  /// Convert JSON exports from PHPMyAdmin to a useable JSON file
// ProcessJSON {
//     /// Input JSON file, exported from PHPMyAdmin
//     #[arg(short, long, verbatim_doc_comment)]
//     in_file: String,
//     /// Key attribute, such as accession, to be used in the new map
//     #[arg(short, long, verbatim_doc_comment)]
//     key: String,
//     /// Optional: Output file. If not provided, will print to stdout
//     #[arg(long, short)]
//     outfile: Option<String>,
// },

// Some(Commands::ProcessJSON {
//     in_file,
//     key,
//     outfile,
// }) => process_json(in_file, key, outfile).expect("JSON Parse Failed"),
