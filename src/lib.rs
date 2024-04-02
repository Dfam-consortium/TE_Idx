use noodles::bgzf;
use serde::{Deserialize, Serialize};
use std::fs::{create_dir_all, File, OpenOptions};
use std::io::{stdout, BufRead, BufReader, Result, Write};
use std::num::NonZeroUsize;
use std::path::Path;
use std::process::exit;
use tempfile::tempfile;

mod idx;

const DATA_DIR: &str = "/home/agray/te_idx/data";
// const ASSEMBLY_DIR: &str = "assembly_alignments";
// const BENCHMARK_DIR: &str = "./data/benchmark_alignments";
// const SEQUENCE_DIR: &str = "./data/sequence";

pub fn bgzf_filter(
    infile: &String,
    position: &usize,
    term: &Option<String>,
    outfile: &Option<String>,
) -> Result<()> {
    if !Path::new(&infile).exists() {
        println!("{} Not Found", &infile);
        std::process::exit(1)
    }

    let worker_count: NonZeroUsize = match NonZeroUsize::new(5) {
        Some(n) => n,
        None => unreachable!(),
    };
    let in_f = File::open(infile).expect("Could Not Open Input File");
    let reader = bgzf::MultithreadedReader::with_worker_count(worker_count, in_f);
    let mut writer: Box<dyn Write> = match outfile {
        Some(outfile) => Box::new(bgzf::MultithreadedWriter::with_worker_count(
            worker_count,
            OpenOptions::new()
                .create(true)
                .append(true)
                .open(&outfile)?,
        )),
        None => Box::new(bgzf::Writer::new(stdout())),
    };

    for result in reader.lines() {
        let line = result?;
        let mut fields: Vec<_> = line.split_whitespace().collect();
        if term.is_none()
            || (fields.len() >= position - 1
                && term.is_some()
                && fields.get(position - 1).unwrap() == term.as_ref().unwrap())
        {
            fields.truncate(14); // TODO readjust this!
            writer
                .write_all(format!("{}\n", fields.join("\t")).as_bytes())
                .expect("Unable to write line");
        }
    }
    Ok(())
}

pub fn prep_beds(in_tsv: &String) -> Result<()> {
    if !Path::new(&in_tsv).exists() {
        eprintln!("Input TSV \"{}\" Not Found", &in_tsv);
        exit(1)
    }

    let chunks: Vec<_> = in_tsv.split("-").collect();
    let assembly = chunks[0];
    let acc_order = chunks[1] == "byacc";
    let bench = chunks[2] == "bench_region.tsv";

    if !acc_order {
        eprintln!("Input TSV \"{}\" Not In Accession Order", &in_tsv);
        exit(1)
    }

    let assemly_dir = format!("{}/{}", &DATA_DIR, &assembly);
    let algin_dir = format!("{}/{}", &assemly_dir, "assembly_alignments");
    let bench_dir = format!("{}/{}", &assemly_dir, "benchmark_alignments");
    let seq_dir = format!("{}/{}", &assemly_dir, "sequences");
    if !Path::new(&assemly_dir).exists() {
        create_dir_all(&assemly_dir)?;
        create_dir_all(&algin_dir)?;
        create_dir_all(&bench_dir)?;
        create_dir_all(&seq_dir)?;
    }

    let target_dir = if bench { bench_dir } else { algin_dir };

    let worker_count: NonZeroUsize = match NonZeroUsize::new(10) {
        Some(n) => n,
        None => unreachable!(),
    };

    let in_f = File::open(in_tsv).expect("Could Not Open Input File");
    let lines = BufReader::new(in_f).lines();
    let mut current_acc = "".to_string();
    let mut out_f = tempfile()?;
    let mut out_writer = bgzf::MultithreadedWriter::with_worker_count(worker_count, out_f);
    for result in lines {
        let line = result?;
        if !&line.starts_with('#') {
            let fields: Vec<_> = line.split_whitespace().collect();
            if fields[3] != current_acc {
                current_acc = fields[3].to_string();
                out_f = File::create(format!("{target_dir}/{current_acc}.bed.bgz",))
                    .expect("Could Not Open Output File");
                out_writer = bgzf::MultithreadedWriter::with_worker_count(worker_count, out_f);
                println!("{}", current_acc);
            };

            out_writer
                .write_all(format!("{line}\n").as_bytes())
                .expect("Unable to write line");
        }
    }

    Ok(())
}

// API Service Subprocesses ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pub fn read_family_assembly_annotations(
    id: &String,
    assembly_id: &String,
    nrph: &bool,
    outfile: &Option<String>,
) {
    let assembly_path: String = format!("{}/{}", &DATA_DIR, &assembly_id);
    if !Path::new(&assembly_path).exists() {
        eprintln!("Assembly \"{}\" Does Not Exist", assembly_path);
        exit(1)
    }
    let fam_file: String = format!("{}/assembly_alignments/{}.bed.bgz", &assembly_path, &id);
    if !Path::new(&fam_file).exists() {
        eprintln!("Family {} Not Found In Assembly {}", id, assembly_path);
        exit(1)
    }

    let position: usize = 13;
    let term: Option<String> = if *nrph { Some("1".to_string()) } else { None };
    match bgzf_filter(&fam_file, &position, &term, outfile) {
        Ok(()) => exit(0),
        Err(err) => {
            eprintln!("Error Filtering File: {} - {}", fam_file, err);
            exit(1)
        }
    }
}

#[derive(Debug)]
#[allow(dead_code)]
#[derive(Serialize, Deserialize)]
struct Annotation {
    accession: String,
    seq_start: String,
    seq_end: String,
    strand: String,
    ali_start: String,
    ali_end: String,
    model_start: String,
    model_end: String,
    bit_score: String,
    e_value: String,
    sequence: String,
    // nrph_hit: String,
    // divergence: String,
    // family_name: String,
    // cigar: String,
    // caf: String,
}
pub fn nhmmer_query(
    assembly: &String,
    chrom: &String,
    start: u64,
    end: u64,
    family: &Option<String>,
    nrph: &bool,
) {
    let assembly_path: String = format!("{}/{}", &DATA_DIR, &assembly);
    // confirm assembly_id and ensure that it accessable
    if !Path::new(&assembly_path).exists() {
        eprintln!("Assembly \"{}\" Does Not Exist", assembly_path);
        exit(1)
    }

    let (filenames, bgz_dir, mut contig_index, index_file) = match idx::prep_idx(&assembly_path) {
        Ok(res) => res,
        Err(e) => panic!("Search Prep Failed, Index may not exist - {:?}", e),
    };

    if !Path::new(&index_file).exists() {
        eprintln!("Assembly \"{}\" Is Not Indexed", assembly_path);
        exit(1)
    }

    let results = idx::search_idx(
        &filenames,
        &bgz_dir,
        &mut contig_index,
        &index_file,
        &chrom,
        start,
        end,
        family,
        *nrph,
        true,
    );

    let mut formatted = Vec::<Annotation>::new();
    match &results {
        Err(e) => {
            eprintln!("Index Search Failed - {}", e);
            exit(1);
        }
        Ok(l) if l.as_slice().is_empty() => {
            println!("No Results Found");
            exit(0);
        }
        Ok(l) => {
            for i in 0..l.len() {
                let fields = l[i].split_whitespace().collect::<Vec<&str>>();
                let annotation = Annotation {
                    accession: fields[3].to_string(),
                    bit_score: fields[4].to_string(),
                    e_value: fields[11].to_string(),
                    model_start: fields[9].to_string(),
                    model_end: fields[10].to_string(),
                    strand: fields[5].to_string(),
                    ali_start: fields[7].to_string(),
                    ali_end: fields[8].to_string(),
                    seq_start: fields[1].to_string(),
                    seq_end: fields[2].to_string(),
                    sequence: fields[0].to_string(),
                };
                formatted.push(annotation)
            }
        }
    };
    match serde_json::to_string(&formatted) {
        Err(e) => {
            eprintln!("Error Converting Results to JSON - {e}");
            exit(1);
        }
        Ok(json_str) => {
            stdout()
                .write_all(json_str.as_bytes())
                .expect("Unable to write line");
            exit(0)
        }
    }
}

#[derive(Serialize, Deserialize)]
struct MaskHit {
    sequence: String,
    seq_start: String,
    seq_end: String,
    repeat_str: String,
    repeat_length: String,
}
pub fn trf_query(assembly: &String, chrom: &String, start: u64, end: u64) -> Result<()> {
    // TODO this method just iterates through all the mask annotations. would be faster to index
    let maskfile = format!("{}/{}/masks/mask.csv.gz", &DATA_DIR, assembly);
    if !Path::new(&maskfile).exists() {
        println!("{} Not Found", &maskfile);
        std::process::exit(1)
    }

    let worker_count: NonZeroUsize = match NonZeroUsize::new(5) {
        Some(n) => n,
        None => unreachable!(),
    };
    let in_f = File::open(maskfile).expect("Could Not Open Input File");
    let reader = bgzf::MultithreadedReader::with_worker_count(worker_count, in_f);
    let mut formatted = Vec::<MaskHit>::new();

    // seq_accession: 0 seq_start: 1 seq_end: 2 repeat_str: 3 repeat_length: 4
    for result in reader.lines() {
        let line = result?;
        let fields: Vec<_> = line.split(",").collect();
        let row_start: u64 = fields[1].parse().expect("Problem with start value");
        let row_end: u64 = fields[2].parse().expect("Problem with end value");
        if &fields[0] == chrom
            && ((row_start >= start && row_start <= end)
                || (row_end <= end && row_end >= start)
                || (row_start < start && row_end > end))
        {
            formatted.push(MaskHit {
                sequence: fields[0].to_string(),
                seq_start: fields[1].to_string(),
                seq_end: fields[2].to_string(),
                repeat_str: fields[3].to_string(),
                repeat_length: fields[4].to_string(),
            })
        }
        match serde_json::to_string(&formatted) {
            Err(e) => {
                eprintln!("Error Converting Results to JSON - {e}");
                exit(1);
            }
            Ok(json_str) => {
                stdout()
                    .write_all(json_str.as_bytes())
                    .expect("Unable to write line");
                exit(0)
            }
        }
    }
    Ok(())
}

pub fn seq_query(assembly: &String, chrom: &String) -> Result<()> {
    let seqfile = format!("{}/{}/sequences/sequence.csv.gz", &DATA_DIR, assembly);
    if !Path::new(&seqfile).exists() {
        println!("{} Not Found", &seqfile);
        std::process::exit(1)
    }

    let worker_count: NonZeroUsize = match NonZeroUsize::new(5) {
        Some(n) => n,
        None => unreachable!(),
    };
    let in_f = File::open(seqfile).expect("Could Not Open Input File");
    let reader = bgzf::MultithreadedReader::with_worker_count(worker_count, in_f);

    let chrom_id: String = "chr".to_owned() + chrom;
    // accession: 0, id: 1, description: 2, length: 3, updated: 4, created: 5, is_genomic: 6
    for result in reader.lines() {
        let line = result?;
        let fields: Vec<_> = line.split(",").collect();
        let seq_id = fields[1];

        if seq_id == chrom || seq_id == &chrom_id {
            stdout()
                .write_all(format!("{}\n", fields[0]).as_bytes())
                .expect("Unable to write line");
        }
    }
    Ok(())
}

pub fn coverage_query(assembly: &String, fam: &String) -> Result<()> {
    let covfile = format!("{}/{}/coverage/coverage_data.tsv.gz", &DATA_DIR, assembly);
    if !Path::new(&covfile).exists() {
        println!("{} Not Found", &covfile);
        std::process::exit(1)
    }

    let worker_count: NonZeroUsize = match NonZeroUsize::new(5) {
        Some(n) => n,
        None => unreachable!(),
    };
    let in_f = File::open(covfile).expect("Could Not Open Input File");
    let reader = bgzf::MultithreadedReader::with_worker_count(worker_count, in_f);
    // family_accession: 0, reversed: 1, forward: 2, nrph: 3, num_full: 4, num_full_nrph: 5, num_rev: 6, karyotype: 7
    for result in reader.lines() {
        let line = result?;
        let fields: Vec<_> = line.split("\t").collect();
        if fields[0] == fam {
            stdout()
            .write_all(format!("{}\n", fields[7]).as_bytes())
            .expect("Unable to write line");
        }
    }
    Ok(())
}

pub fn percent_id_query(){
    // family_accession,threshold,graph_json,max_insert,num_seqs
    todo!()
}