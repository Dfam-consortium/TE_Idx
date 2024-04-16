use noodles::bgzf;
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use std::collections::HashMap;
use std::fs::{create_dir_all, read_to_string, File, OpenOptions};
use std::io::{stdout, BufRead, BufReader, Result, Write};
use std::num::NonZeroUsize;
use std::path::Path;
use std::process::exit;
use tempfile::tempfile;

mod idx;

pub const DATA_DIR: &'static str = "/home/agray/te_idx/data";
pub const ASSEMBLY_DIR: &'static str = "assembly_alignments";
pub const BENCHMARK_DIR: &'static str = "benchmark_alignments";
pub const MASKS_DIR: &'static str = "masks";
pub const MOD_LEN_DIR: &'static str = "model_lengths";
pub const SEQUENCE_DIR: &'static str = "sequences";

// .bed fields => seq_id 0, seq_start 1, seq_end 2, family_accession 3, hit_bit_score 4, strand 5, ali_start 6, ali_end 7,
//                model_start 8, model_end 9, hit_evalue_score 10, nrph_hit 11, divergence 12, family_name 13, cigar 14, caf 15

// fam-annotations => sequence name 0, model accession 1, model name 2, bit score 3, e-value 4, hmm start 5, hmm end 6,
//                   hmm length 7, strand 8, alignment start 9, alignment end 10, envelope start 11, envelope end 12, sequence length 13

// .bed to fam-annotations:
// 0->0, from seq file
// 3->1
// 13->2
// 4->3
// 10->4
// 8->5
// 9->6
// 3->7 from family table
// 5->8
// 1->9
// 2->10
// 6->11
// 7->12
// 0->13 from seq file

trait Formattable {
    fn to_json(&self) -> serde_json::Value;
}

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

#[derive(Serialize, Deserialize)]
struct MaskHit {
    sequence: String,
    seq_start: String,
    seq_end: String,
    repeat_str: String,
    repeat_length: String,
}

impl Formattable for MaskHit {
    fn to_json(&self) -> serde_json::Value {
        json!({
            "sequence": self.sequence,
            "seq_start": self.seq_start,
            "seq_end": self.seq_end,
            "repeat_str": self.repeat_str,
            "repeat_length": self.repeat_length,
        })
    }
}

impl Formattable for Annotation {
    fn to_json(&self) -> serde_json::Value {
        json!({
            "accession": self.accession,
            "seq_start": self.seq_start,
            "seq_end": self.seq_end,
            "strand": self.strand,
            "ali_start": self.ali_start,
            "ali_end": self.ali_end,
            "model_start": self.model_start,
            "model_end": self.model_end,
            "bit_score": self.bit_score,
            "e_value": self.e_value,
            "sequence": self.sequence,
        })
    }
}

fn build_annotation(fields: Vec<&str>) -> Box<dyn Formattable> {
    Box::new(Annotation {
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
    })
}

fn build_mask(fields: Vec<&str>) -> Box<dyn Formattable> {
    Box::new(MaskHit {
        sequence: fields[0].to_string(),
        seq_start: fields[1].to_string(),
        seq_end: fields[2].to_string(),
        repeat_str: fields[3].to_string(),
        repeat_length: fields[4].to_string(),
    })
}

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

pub fn prep_beds(assembly: &String, in_tsv: &String, data_type: &String) -> Result<()> {
    if !Path::new(&in_tsv).exists() {
        eprintln!("Input TSV \"{}\" Not Found", &in_tsv);
        exit(1)
    }

    let db_dir = format!("{}/{}", &DATA_DIR, &assembly);
    let target_dir = format!("{}/{}", &db_dir, &data_type);
    if !Path::new(&db_dir).exists() {
        create_dir_all(&db_dir)?;
    }

    if !Path::new(&target_dir).exists() {
        create_dir_all(&target_dir)?;
    }

    let acc_idx: usize = match data_type.as_str() {
        ASSEMBLY_DIR => 3,
        BENCHMARK_DIR => 3,
        MASKS_DIR => 0,
        SEQUENCE_DIR => 0,
        _ => panic!("Invalid Data Type"),
    };

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
            if fields[acc_idx] != current_acc {
                // assume accession order TODO confirm this
                current_acc = fields[acc_idx].to_string();
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

pub fn process_json(in_file: &String, key: &String, outfile: &Option<String>) -> Result<()> {
    if !Path::new(&in_file).exists() {
        println!("{} Not Found", &in_file);
        std::process::exit(1)
    }
    let in_f = File::open(in_file).expect("Could Not Open Input File");

    let in_data: Value = serde_json::from_reader(in_f)?;

    let mut out_data = HashMap::new();
    if let Some(items) = in_data[2]["data"].as_array() {
        for item in items {
            let mut new_item = item.clone().as_object().unwrap().to_owned();
            new_item.remove(key);
            out_data.insert(item[key].to_string().replace("\"", ""), new_item);
        }
    }
    let json = json!({
        "assembly": in_data[1]["name"],
        "version": in_data[0]["version"],
        "data": out_data
    });

    let output = serde_json::to_string(&json)?;
    let mut writer: Box<dyn Write> = match outfile {
        Some(outfile) => Box::new(File::create(&outfile).expect("Could Not Open Output File")),
        None => Box::new(stdout()),
    };
    writer
        .write_all(format!("{}\n", output).as_bytes())
        .expect("Unable to write JSON");
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
    let fam_file: String = format!("{}/{}/{}.bed.bgz", &assembly_path, &ASSEMBLY_DIR, &id);
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

pub fn idx_query(
    assembly: &String,
    data_type: &String,
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

    let (filenames, bgz_dir, mut contig_index, index_file) =
        match idx::prep_idx(&assembly_path, &data_type) {
            Ok(res) => res,
            Err(e) => panic!("Search Prep Failed, Index may not exist - {:?}", e),
        };

    if !Path::new(&index_file).exists() {
        eprintln!(
            "Assembly \"{}\" Is Not Indexed For {}",
            assembly_path, &data_type
        );
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

    let formatter = match data_type.as_str() {
        ASSEMBLY_DIR => build_annotation,
        MASKS_DIR => build_mask,
        _ => panic!("Invalid formatter type"),
    };
    let mut formatted = Vec::new();
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
                formatted.push(formatter(fields).to_json());
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

pub fn json_query(
    assembly: &String,
    data_type: &String,
    key: &String,
    target: &String,
) -> Result<String> {
    let target_file = format!(
        "{}/{}/{}/{}-{}-processed.json",
        &DATA_DIR, &assembly, &data_type, &assembly, &data_type
    );
    if !Path::new(&target_file).exists() {
        println!("{} Not Found", &target_file);
        std::process::exit(1)
    }

    let in_str = read_to_string(&target_file).expect("Could Not Read String");
    let in_data: Value = serde_json::from_str(&in_str).expect("JSON was not well-formatted");

    let val = in_data
        .get("data")
        .and_then(|data| data.get(key).and_then(|item| item.get(target)))
        .expect("Key Target Pair Not Found");
    let ret = val.as_str().expect("asd");
    return Ok(ret.to_string());
}
