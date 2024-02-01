use noodles::bgzf;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fs::File;
use std::io::{stdout, BufRead, Result, Write};
use std::num::NonZeroUsize;
use std::path::Path;
use std::process::exit;

mod idx;

const DATA_DIR: &str = "./data";
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
            File::create(outfile).expect("Could Not Open Output File"),
        )),
        None => Box::new(bgzf::Writer::new(stdout())),
    };

    for line in reader.lines() {
        // TODO replace unwrap with match
        let linestr = line.unwrap();
        let mut fields: Vec<_> = linestr.split_whitespace().collect();
        if term.is_none() || fields.get(position - 1).unwrap() == term.as_ref().unwrap() {
            fields.truncate(9);
            writer
                .write_all(format!("{}\n", fields.join("\t")).as_bytes())
                .expect("Unable to write line");
        }
    }
    Ok(())
}

pub fn read_family_assembly_annotations(id: &String, assembly_id: &String, nrph: &bool) {
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
    match bgzf_filter(&fam_file, &position, &term, &None) {
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
    family_accession: String,
    seq_start: String,
    seq_end: String,
    strand: String,
    // ali_start: String,
    // ali_end: String,
    model_start: String,
    model_end: String,
    hit_bit_score: String,
    hit_evalue_score: String,
    nrph_hit: String,
    chrom: String,
}
pub fn read_annotations(
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
    );

    let mut field_map = HashMap::new();
    field_map.insert("family_accession", 3);
    field_map.insert("seq_start", 1);
    field_map.insert("seq_end", 2);
    field_map.insert("strand", 5);
    // field_map.insert("ali_start", 4);
    // field_map.insert("ali_end", 4);
    field_map.insert("model_start", 7);
    field_map.insert("model_end", 8);
    field_map.insert("hit_bit_score", 4);
    field_map.insert("hit_evalue_score", 6);
    field_map.insert("nrph_hit", 12);
    field_map.insert("chrom", 0);

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
                    family_accession: fields[*field_map.get("family_accession").unwrap()]
                        .to_string(),
                    seq_start: fields[*field_map.get("seq_start").unwrap()].to_string(),
                    seq_end: fields[*field_map.get("seq_end").unwrap()].to_string(),
                    strand: fields[*field_map.get("strand").unwrap()].to_string(),
                    // ali_start: fields[*field_map.get("ali_start").unwrap()].to_string(),
                    // ali_end: fields[*field_map.get("ali_end").unwrap()].to_string(),
                    model_start: fields[*field_map.get("model_start").unwrap()].to_string(),
                    model_end: fields[*field_map.get("model_end").unwrap()].to_string(),
                    hit_bit_score: fields[*field_map.get("hit_bit_score").unwrap()].to_string(),
                    hit_evalue_score: fields[*field_map.get("hit_evalue_score").unwrap()]
                        .to_string(),
                    nrph_hit: fields[*field_map.get("nrph_hit").unwrap()].to_string(),
                    chrom: fields[*field_map.get("chrom").unwrap()].to_string(),
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

// pub fn find_family(_id: &String, assembly: &String) {
//     let paths = read_dir(format!("{}/{}/{}", DATA_DIR, assembly, ASSEMBLY_DIR)).unwrap();
//     for path in paths {
//         // TODO
//         // if path.ends_with(format!("{}.bed.bgz", id)) {
//         println!("Name: {}", path.unwrap().path().display())
//         // }
//     }
// }
