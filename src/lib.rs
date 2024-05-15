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

pub const DATA_DIR: &'static str = "/home/agray/te_idx/data"; // TODO replace with config
pub const EXPORT_DIR: &'static str = "/home/agray/te_idx/exports"; // TODO replace with config

pub const ASSEMBLY_DIR: &'static str = "assembly_alignments";
pub const ASSEMBLY_FILE: &'static str = "-byacc-full_region.tsv";
pub const BENCHMARK_DIR: &'static str = "benchmark_alignments";
pub const BENCHMARK_FILE: &'static str = "-byacc-bench_region.tsv";
pub const MASKS_DIR: &'static str = "masks";
pub const MASKS_FILE: &'static str = "-mask.tsv";
pub const MOD_LEN_DIR: &'static str = "model_lengths";
pub const MOD_LEN_FILE: &'static str = "-model_lengths.json";
pub const SEQUENCE_DIR: &'static str = "sequences";
pub const SEQUENCE_FILE: &'static str = "-processed-sequences.json";

const DATA_ELEMENTS: [&str; 5] = [
    ASSEMBLY_DIR,
    BENCHMARK_DIR,
    MASKS_DIR,
    MOD_LEN_DIR,
    SEQUENCE_DIR,
];
pub const INDEX_DATA_TYPES: [&str; 3] = [ASSEMBLY_DIR, BENCHMARK_DIR, MASKS_DIR];
pub const JSON_DATA_TYPES: [&str; 2] = [MOD_LEN_DIR, SEQUENCE_DIR];

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

fn download_format<'a>(
    fields: &Vec<&'a str>,
    seq_name: &'a str,
    hmm_len: &'a str,
    seq_len: &'a str,
) -> Vec<&'a str> {
    let mut fmtted: Vec<&'a str> = Vec::with_capacity(fields.len());
    fmtted.push(seq_name); // seq name
    fmtted.push(fields[3]); // Family Accession
    fmtted.push(fields[14]); // Family name
    fmtted.push(fields[4]); // bit score
    fmtted.push(fields[11]); // e-value
    fmtted.push(fields[9]); // hmm start
    fmtted.push(fields[10]); // hmm end
    fmtted.push(hmm_len); // hmm length
    fmtted.push(fields[5]); // strand
    fmtted.push(fields[7]); // align start
    fmtted.push(fields[8]); // align end
    fmtted.push(fields[1]); // env start
    fmtted.push(fields[2]); // env end
    fmtted.push(seq_len); // sequence length
    return fmtted;
}

pub fn bgzf_filter(
    assembly: &String,
    data_type: &String,
    fam: &String,
    position: &usize,
    term: &Option<String>,
    outfile: &Option<String>,
    dl_fmt: bool,
    data_directory: Option<&str>,
) -> Result<()> {
    let data_dir = data_directory.unwrap_or(DATA_DIR);
    let assembly_path: String = format!("{}/{}/{}", &data_dir, &assembly, &data_type);
    if !Path::new(&assembly_path).exists() {
        eprintln!("Data \"{}\" Does Not Exist", assembly_path);
        exit(1)
    }
    let fam_file: String = format!("{}/{}.bed.bgz", &assembly_path, &fam);
    if !Path::new(&fam_file).exists() {
        eprintln!("Family {} Not Found In Assembly {}", &fam, assembly_path);
        exit(1)
    }

    let worker_count: NonZeroUsize = match NonZeroUsize::new(5) {
        Some(n) => n,
        None => unreachable!(),
    };
    let in_f = File::open(fam_file).expect("Could Not Open Input File");
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

    let header;
    if dl_fmt {
        header = "#sequence name	model accession	model name	bit score	e-value	hmm start	hmm end	hmm length	strand	alignment start	alignment end	envelope start	envelope end	sequence length";
    } else {
        header = "#seq_id\tseq_start\tseq_end\tfamily_accession\thit_bit_score\tstrand\tbias\tali_start\tali_end\tmodel_start\tmodel_end\thit_evalue_score\tnrph_hit\tdivergence\t*family_name\t*cigar\t*caf";
    }
    writer
        .write_all(format!("{}\n", header).as_bytes())
        .expect("Unable to write line");

    let hmm_len = match json_query(
        &assembly,
        &MOD_LEN_DIR.to_string(),
        &fam.to_string(),
        &"length".to_string(),
        None,
    ) {
        Ok(len) => len,
        Err(e) => panic!("{}", e),
    };

    let in_str = read_to_string(&format!(
        "{}/{}/{}/{}{}",
        &data_dir, &assembly, &SEQUENCE_DIR, &assembly, &SEQUENCE_FILE
    ))
    .expect("Could Not Read String");
    let seq_json: Value = serde_json::from_str(&in_str).expect("JSON was not well-formatted");
    let seq_data = seq_json
        .get("data")
        .expect(&format!("Sequence Info For {} Not Found", &fam));

    let mut seq_name;
    let mut seq_len;

    for result in reader.lines() {
        let line = result?;
        let mut fields: Vec<_> = line.split_whitespace().collect();
        if term.is_none()
            || (fields.len() >= position - 1
                && term.is_some()
                && fields.get(position - 1).unwrap() == term.as_ref().unwrap())
        {
            if dl_fmt {
                let chrom_id = &fields[0].to_string();
                let seq_details = seq_data
                    .get(chrom_id)
                    .expect(&format!("Sequence {} Not Found", chrom_id));
                seq_name = seq_details
                    .get(&"id".to_string())
                    .expect(&format!("Name Not Found For {}", chrom_id))
                    .as_str()
                    .expect("Couldn't Cast To Str");
                seq_len = seq_details
                    .get(&"length".to_string())
                    .expect(&format!("Length Not Found For {}", chrom_id))
                    .as_str()
                    .expect("Couldn't Cast To Str");
                fields = download_format(&fields, seq_name, &hmm_len, seq_len);
            } else {
                fields.truncate(14); // TODO readjust this!
            }
            writer
                .write_all(format!("{}\n", &fields.join("\t")).as_bytes())
                .expect("Unable to write line");
        }
    }
    Ok(())
}

// Setup Methods ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pub fn prep_beds(
    assembly: &String,
    in_tsv: &String,
    data_type: &String,
    data_directory: Option<&str>,
) -> Result<()> {
    let data_dir = data_directory.unwrap_or(DATA_DIR);
    if !Path::new(&in_tsv).exists() {
        eprintln!("Input TSV \"{}\" Not Found", &in_tsv);
        exit(1)
    }

    let db_dir = format!("{}/{}", &data_dir, &assembly);
    let target_dir = format!("{}/{}", &db_dir, &data_type);
    if !Path::new(&db_dir).exists() {
        create_dir_all(&db_dir)?;
    }

    if !Path::new(&target_dir).exists() {
        create_dir_all(&target_dir)?;
    }

    let acc_idx: usize = match data_type.as_str() {
        ASSEMBLY_DIR => 3,
        BENCHMARK_DIR => 1,
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
                // println!("{}", current_acc);
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
        eprintln!("{} Not Found", &in_file);
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

pub fn prepare_assembly(
    assembly: &String,
    data_directory: Option<&str>,
    export_directory: Option<&str>,
) -> Result<()> {
    let data_dir = data_directory.unwrap_or(DATA_DIR);
    let export_dir = export_directory.unwrap_or(EXPORT_DIR);
    if !Path::new(&data_dir).exists() {
        eprintln!("{} Not Found", &data_dir);
        std::process::exit(1)
    }
    if !Path::new(&export_dir).exists() {
        eprintln!("{} Not Found", &export_dir);
        std::process::exit(1)
    }

    let export_dir = format!("{}/{}", &export_dir, &assembly);
    let working_dir = format!("{}/{}", &data_dir, &assembly);
    if !Path::new(&export_dir).exists() {
        eprintln!("Assembly Export Not Found - {}", &export_dir);
        std::process::exit(1)
    }
    if !Path::new(&working_dir).exists() {
        println!(
            "Target Assembly Directory Not Found, Creating {},",
            &working_dir
        );
        create_dir_all(&working_dir)?;
    }

    fn file_to_source(s: &str) -> Option<&'static str> {
        match s {
            ASSEMBLY_DIR => Some(ASSEMBLY_FILE),
            BENCHMARK_DIR => Some(BENCHMARK_FILE),
            MASKS_DIR => Some(MASKS_FILE),
            MOD_LEN_DIR => Some(MOD_LEN_FILE),
            SEQUENCE_DIR => Some(SEQUENCE_FILE),
            _ => None,
        }
    }
    let mut planner = HashMap::new();
    for element in DATA_ELEMENTS {
        let source = format!(
            "{}/{}{}",
            export_dir,
            &assembly,
            file_to_source(&element).unwrap().to_string()
        );
        let target = format!("{}/{}", &working_dir, &element);

        let have_source = Path::new(&source).exists();
        let have_target =
            Path::new(&target).exists() && !Path::new(&target).read_dir()?.next().is_none();
        let needed = have_source && !have_target;

        let info = HashMap::from([
            ("source", source),
            ("target", target),
            ("needed", needed.to_string()),
        ]);
        planner.insert(element, info);
    }

    for element in DATA_ELEMENTS {
        if matches!(
            planner
                .get(element)
                .unwrap()
                .get("needed")
                .unwrap()
                .as_str(),
            "true"
        ) {
            println!("Preparing {}: ", element);
            let target = planner.get(element).unwrap().get("target").unwrap();
            if !Path::new(target).exists() {
                println!("   Target Directory Not Found, Creating {},", &target);
                create_dir_all(&target)?;
            }
            let source = planner.get(element).unwrap().get("source").unwrap();
            if source.ends_with(".json") {
                let key = match element {
                    MOD_LEN_DIR => "family_accession",
                    SEQUENCE_DIR => "accession",
                    _ => "accession",
                };
                let outfile = format!(
                    "{}/{}-{}{}",
                    target,
                    assembly,
                    "processed",
                    file_to_source(element).unwrap().to_string()
                );
                process_json(source, &key.to_string(), &Some(outfile))
                    .expect(&format!("Processing {} JSON Failed", element));
            } else if source.ends_with(".tsv") {
                println!("   Splitting And Compressing BED Files For {}", element);
                prep_beds(assembly, source, &element.to_string(), None)
                    .expect("BED File Prep Failed");
                println!("   Indexing {}", element);
                let (filenames, bgz_dir, mut contig_index, index_file) =
                    idx::prep_idx(&working_dir, &element.to_string()).expect("Index Prep Failed");
                idx::build_idx(&filenames, &bgz_dir, &mut contig_index, &index_file)
                    .expect("Indexing Failed")
            } else {
                eprintln!("Source type not recognized - {}", source)
            }
        }
    }
    Ok(())
}

// API Service Subprocesses ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pub fn read_family_assembly_annotations(
    id: &String,
    assembly_id: &String,
    nrph: &bool,
    outfile: &Option<String>,
    data_directory: Option<&str>,
) {
    let data_dir = data_directory.unwrap_or(DATA_DIR);
    let assembly_path: String = format!("{}/{}", &data_dir, &assembly_id);
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
    match bgzf_filter(
        &assembly_id,
        &ASSEMBLY_DIR.to_string(),
        &id,
        &position,
        &term,
        outfile,
        true,
        None,
    ) {
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
    data_directory: Option<&str>,
) -> Result<()> {
    let data_dir = data_directory.unwrap_or(DATA_DIR);
    let assembly_path: String = format!("{}/{}", &data_dir, &assembly);
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
    data_directory: Option<&str>,
) -> Result<String> {
    let data_dir = data_directory.unwrap_or(DATA_DIR);
    let target_file = format!(
        "{}/{}/{}/{}-processed-{}.json",
        &data_dir, &assembly, &data_type, &assembly, &data_type
    );
    if !Path::new(&target_file).exists() {
        eprintln!("{} Not Found", &target_file);
        std::process::exit(1)
    }

    let in_str = read_to_string(&target_file).expect("Could Not Read String");
    let in_data: Value = serde_json::from_str(&in_str).expect("JSON was not well-formatted");

    let val = in_data
        .get("data")
        .and_then(|data| data.get(key).and_then(|item| item.get(target)))
        .expect("Key Target Pair Not Found");
    let ret = val.as_str().expect("Can't String to str");
    return Ok(ret.to_string());
}
