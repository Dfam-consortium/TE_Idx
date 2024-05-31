use noodles::bgzf;
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use std::collections::HashMap;
use std::fs::{copy, create_dir_all, read_to_string, File, OpenOptions};
use std::io::{stdout, BufRead, BufReader, Result, Write};
use std::num::NonZeroUsize;
use std::path::Path;
use tempfile::tempfile;

pub mod idx;

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
pub const SEQUENCE_FILE: &'static str = "-sequences.json";

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
    fn from_export_tsv(tsv_line: &Vec<&str>) -> Self;
    fn from_bed(bed_line: &Vec<&str>) -> Self;
    fn to_json(&self) -> serde_json::Value;
    fn to_bed_fmt(&self) -> Vec<&str>;
    fn to_dl_fmt(&self, seq_name: &str, hmm_len: &str) -> Vec<String>;
    fn to_filter_fmt(&self) -> Vec<String>;
    fn get_acc(&self) -> String;
}

#[derive(Serialize, Deserialize)]
pub struct Annotation {
    seq_acc: String,     // Dfamseq accession for sequences in assembly (1..)
    fam_acc: String,     // Dfam family accession ( e.g DF######### )
    family_name: String, // Optional family name
    bit_score: String,   // Alignment bitscore
    e_value: String,     // Alignment evalue
    bias: String,        // God knows...(actually it's in the nhmmer manual)
    model_start: String, // pHMM start position (1-based, fully closed)
    model_end: String,   // pHMM end position (1-based, fully closed)
    strand: String,      // '+' or '-'
    ali_start: String,   // nhmmer "envelope" start
    ali_end: String,     // nhmmer "envelope" end
    seq_start: String,   // Dfamseq sequence start
    seq_end: String,     // Dfamseq sequence end
    seq_len: String,     // The length of the Dfamseq
    cigar: String,       // CIGAR string of sequence alignment
    kimura_div: String, // Kimura percent divergence ( only in full_region, not in benchmark_region )
    nrph_hit: String,   // NRPH (1 or 0) (only in full_region, not in benchmark_region)
    caf: String,        // Compressed Alignment Format (CAF) (only in full_region...)
}
impl Formattable for Annotation {
    fn from_export_tsv(tsv_line: &Vec<&str>) -> Self {
        Self {
            seq_acc: tsv_line[0].to_string(),
            fam_acc: tsv_line[1].to_string(),
            family_name: tsv_line[2].to_string(),
            bit_score: tsv_line[3].to_string(),
            e_value: tsv_line[4].to_string(),
            bias: tsv_line[5].to_string(),
            model_start: tsv_line[6].to_string(),
            model_end: tsv_line[7].to_string(),
            strand: tsv_line[8].to_string(),
            ali_start: tsv_line[9].to_string(),
            ali_end: tsv_line[10].to_string(),
            seq_start: tsv_line[11].to_string(),
            seq_end: tsv_line[12].to_string(),
            seq_len: tsv_line[13].to_string(),
            cigar: tsv_line[14].to_string(),
            kimura_div: tsv_line[15].to_string(),
            nrph_hit: tsv_line[16].to_string(),
            caf: tsv_line[17].to_string(),
        }
    }

    fn to_bed_fmt(&self) -> Vec<&str> {
        vec![
            &self.seq_acc,
            &self.seq_start,
            &self.seq_end,
            &self.fam_acc,
            &self.bit_score,
            &self.strand,
            &self.bias,
            &self.ali_start,
            &self.ali_end,
            &self.model_start,
            &self.model_end,
            &self.e_value,
            &self.nrph_hit,
            &self.kimura_div,
            &self.family_name,
            &self.seq_len,
            &self.cigar,
            &self.caf,
        ]
    }

    fn from_bed(bed_line: &Vec<&str>) -> Self {
        Self {
            seq_acc: bed_line[0].to_string(),
            seq_start: bed_line[1].to_string(),
            seq_end: bed_line[2].to_string(),
            fam_acc: bed_line[3].to_string(),
            bit_score: bed_line[4].to_string(),
            strand: bed_line[5].to_string(),
            bias: bed_line[6].to_string(),
            ali_start: bed_line[7].to_string(),
            ali_end: bed_line[8].to_string(),
            model_start: bed_line[9].to_string(),
            model_end: bed_line[10].to_string(),
            e_value: bed_line[11].to_string(),
            nrph_hit: bed_line[12].to_string(),
            kimura_div: bed_line[13].to_string(),
            family_name: bed_line[14].to_string(),
            seq_len: bed_line[15].to_string(),
            cigar: bed_line[16].to_string(),
            caf: bed_line[17].to_string(),
        }
    }

    fn to_json(&self) -> serde_json::Value {
        json!({
            "sequence": self.seq_acc,
            "accession": self.fam_acc,
            "bit_score": self.bit_score,
            "e_value": self.e_value,
            "seq_start": self.seq_start,
            "seq_end": self.seq_end,
            "strand": self.strand,
            "ali_start": self.ali_start,
            "ali_end": self.ali_end,
            "model_start": self.model_start,
            "model_end": self.model_end,
            "bit_score": self.bit_score,
            "e_value": self.e_value,
        })
    }

    fn to_dl_fmt(&self, seq_name: &str, hmm_len: &str) -> Vec<String> {
        vec![
            seq_name.to_string(),
            self.fam_acc.clone(),
            self.family_name.clone(),
            self.bit_score.clone(),
            self.e_value.clone(),
            self.model_start.clone(),
            self.model_end.clone(),
            hmm_len.to_string(),
            self.strand.clone(),
            self.ali_start.clone(),
            self.ali_end.clone(),
            self.seq_start.clone(),
            self.seq_end.clone(),
            self.seq_len.clone(),
        ]
    }

    fn to_filter_fmt(&self) -> Vec<String> {
        vec![
            self.seq_acc.to_string(),
            self.seq_start.to_string(),
            self.seq_end.to_string(),
            self.fam_acc.to_string(),
            self.bit_score.to_string(),
            self.strand.to_string(),
            self.bias.to_string(),
            self.ali_start.to_string(),
            self.ali_end.to_string(),
            self.model_start.to_string(),
            self.model_end.to_string(),
            self.e_value.to_string(),
            self.nrph_hit.to_string(),
            self.kimura_div.to_string(),
            self.family_name.to_string(),
            self.seq_len.to_string(),
        ]
    }

    fn get_acc(&self) -> String {
        self.fam_acc.clone()
    }
}

#[derive(Serialize, Deserialize)]
pub struct BenchMarkAnnotation {
    seq_acc: String,     // Dfamseq accession for sequences in assembly (1..)
    fam_acc: String,     // Dfam family accession ( e.g DF######### )
    family_name: String, // Optional family name
    bit_score: String,   // Alignment bitscore
    e_value: String,     // Alignment evalue
    bias: String,        // God knows...(actually it's in the nhmmer manual)
    model_start: String, // pHMM start position (1-based, fully closed)
    model_end: String,   // pHMM end position (1-based, fully closed)
    strand: String,      // '+' or '-'
    ali_start: String,   // nhmmer "envelope" start
    ali_end: String,     // nhmmer "envelope" end
    seq_start: String,   // Dfamseq sequence start
    seq_end: String,     // Dfamseq sequence end
    seq_len: String,     // The length of the Dfamseq
    cigar: String,       // CIGAR string of sequence alignment
}

impl Formattable for BenchMarkAnnotation {
    fn from_export_tsv(tsv_line: &Vec<&str>) -> Self {
        Self {
            seq_acc: tsv_line[0].to_string(),
            fam_acc: tsv_line[1].to_string(),
            family_name: tsv_line[2].to_string(),
            bit_score: tsv_line[3].to_string(),
            e_value: tsv_line[4].to_string(),
            bias: tsv_line[5].to_string(),
            model_start: tsv_line[6].to_string(),
            model_end: tsv_line[7].to_string(),
            strand: tsv_line[8].to_string(),
            ali_start: tsv_line[9].to_string(),
            ali_end: tsv_line[10].to_string(),
            seq_start: tsv_line[11].to_string(),
            seq_end: tsv_line[12].to_string(),
            seq_len: tsv_line[13].to_string(),
            cigar: tsv_line[14].to_string(),
        }
    }

    fn to_json(&self) -> serde_json::Value {
        json!({
            "sequence": self.seq_acc,
            "accession": self.fam_acc,
            "bit_score": self.bit_score,
            "e_value": self.e_value,
            "seq_start": self.seq_start,
            "seq_end": self.seq_end,
            "strand": self.strand,
            "ali_start": self.ali_start,
            "ali_end": self.ali_end,
            "model_start": self.model_start,
            "model_end": self.model_end,
            "bit_score": self.bit_score,
            "e_value": self.e_value,
        })
    }

    fn to_bed_fmt(&self) -> Vec<&str> {
        vec![
            &self.seq_acc,
            &self.seq_start,
            &self.seq_end,
            &self.fam_acc,
            &self.bit_score,
            &self.strand,
            &self.bias,
            &self.ali_start,
            &self.ali_end,
            &self.model_start,
            &self.model_end,
            &self.e_value,
            &self.family_name,
            &self.seq_len,
            &self.cigar,
        ]
    }

    fn from_bed(bed_line: &Vec<&str>) -> Self {
        Self {
            seq_acc: bed_line[0].to_string(),
            seq_start: bed_line[1].to_string(),
            seq_end: bed_line[2].to_string(),
            fam_acc: bed_line[3].to_string(),
            bit_score: bed_line[4].to_string(),
            strand: bed_line[5].to_string(),
            bias: bed_line[6].to_string(),
            ali_start: bed_line[7].to_string(),
            ali_end: bed_line[8].to_string(),
            model_start: bed_line[9].to_string(),
            model_end: bed_line[10].to_string(),
            e_value: bed_line[11].to_string(),
            family_name: bed_line[13].to_string(),
            seq_len: bed_line[14].to_string(),
            cigar: bed_line[15].to_string(),
        }
    }

    fn to_dl_fmt(&self, seq_name: &str, hmm_len: &str) -> Vec<String> {
        vec![
            seq_name.to_string(),
            self.fam_acc.clone(),
            self.family_name.clone(),
            self.bit_score.clone(),
            self.e_value.clone(),
            self.model_start.clone(),
            self.model_end.clone(),
            hmm_len.to_string(),
            self.strand.clone(),
            self.ali_start.clone(),
            self.ali_end.clone(),
            self.seq_start.clone(),
            self.seq_end.clone(),
            self.seq_len.clone(),
        ]
    }

    fn to_filter_fmt(&self) -> Vec<String> {
        vec![
            self.seq_acc.to_string(),
            self.seq_start.to_string(),
            self.seq_end.to_string(),
            self.fam_acc.to_string(),
            self.bit_score.to_string(),
            self.strand.to_string(),
            self.bias.to_string(),
            self.ali_start.to_string(),
            self.ali_end.to_string(),
            self.model_start.to_string(),
            self.model_end.to_string(),
            self.e_value.to_string(),
            self.family_name.to_string(),
            self.seq_len.to_string(),
        ]
    }

    fn get_acc(&self) -> String {
        self.fam_acc.clone()
    }
}

#[derive(Serialize, Deserialize)]
struct MaskHit {
    seq_acc: String,
    seq_start: String,
    seq_end: String,
    repeat_str: String,
    repeat_length: String,
}

impl Formattable for MaskHit {
    fn from_export_tsv(tsv_line: &Vec<&str>) -> Self {
        Self {
            seq_acc: tsv_line[0].to_string(),
            seq_start: tsv_line[1].to_string(),
            seq_end: tsv_line[2].to_string(),
            repeat_str: tsv_line[3].to_string(),
            repeat_length: tsv_line[4].to_string(),
        }
    }
    fn to_json(&self) -> serde_json::Value {
        json!({
            "seq_acc": self.seq_acc,
            "seq_start": self.seq_start,
            "seq_end": self.seq_end,
            "repeat_str": self.repeat_str,
            "repeat_length": self.repeat_length,
        })
    }

    fn to_bed_fmt(&self) -> Vec<&str> {
        vec![
            &self.seq_acc,
            &self.seq_start,
            &self.seq_end,
            &self.repeat_str,
            &self.repeat_length,
        ]
    }

    fn from_bed(bed_line: &Vec<&str>) -> Self {
        MaskHit::from_export_tsv(&bed_line)
    }

    fn to_dl_fmt(&self, _seq_name: &str, _hmm_len: &str) -> Vec<String> {
        vec![
            self.seq_acc.clone(),
            self.seq_start.clone(),
            self.seq_end.clone(),
            self.repeat_str.clone(),
            self.repeat_length.clone(),
        ]
    }

    fn to_filter_fmt(&self) -> Vec<String> {
        vec![
            self.seq_acc.to_string(),
            self.seq_start.to_string(),
            self.seq_end.to_string(),
            self.repeat_str.to_string(),
            self.repeat_length.to_string(),
        ]
    }

    fn get_acc(&self) -> String {
        self.seq_acc.clone()
    }
}

enum FormattableLine {
    Annotation(Annotation),
    BenchMarkAnnotation(BenchMarkAnnotation),
    MaskHit(MaskHit),
}
impl FormattableLine {
    fn from_export_tsv(tsv_line: &Vec<&str>, data_type: &str) -> Self {
        match data_type {
            ASSEMBLY_DIR => FormattableLine::Annotation(Annotation::from_export_tsv(tsv_line)),
            BENCHMARK_DIR => {
                FormattableLine::BenchMarkAnnotation(BenchMarkAnnotation::from_export_tsv(tsv_line))
            }
            MASKS_DIR => FormattableLine::MaskHit(MaskHit::from_export_tsv(tsv_line)),
            _ => panic!("Can't Format!"),
        }
    }

    fn from_bed(bed_line: &Vec<&str>, data_type: &str) -> Self {
        match data_type {
            ASSEMBLY_DIR => FormattableLine::Annotation(Annotation::from_bed(bed_line)),
            BENCHMARK_DIR => {
                FormattableLine::BenchMarkAnnotation(BenchMarkAnnotation::from_bed(bed_line))
            }
            MASKS_DIR => FormattableLine::MaskHit(MaskHit::from_bed(bed_line)),
            _ => panic!("Can't Format!"),
        }
    }

    fn to_dl_fmt(&self, seq_name: &str, hmm_len: &str) -> Vec<String> {
        match self {
            FormattableLine::Annotation(annotation) => annotation.to_dl_fmt(seq_name, hmm_len),
            FormattableLine::BenchMarkAnnotation(benchmark) => {
                benchmark.to_dl_fmt(seq_name, hmm_len)
            }
            FormattableLine::MaskHit(mask_hit) => mask_hit.to_dl_fmt(seq_name, hmm_len),
        }
    }

    fn to_json(&self) -> serde_json::Value {
        match self {
            FormattableLine::Annotation(annotation) => annotation.to_json(),
            FormattableLine::BenchMarkAnnotation(benchmark) => benchmark.to_json(),
            FormattableLine::MaskHit(mask_hit) => mask_hit.to_json(),
        }
    }

    fn to_bed_fmt(&self) -> Vec<&str> {
        match self {
            FormattableLine::Annotation(annotation) => annotation.to_bed_fmt(),
            FormattableLine::BenchMarkAnnotation(benchmark) => benchmark.to_bed_fmt(),
            FormattableLine::MaskHit(mask_hit) => mask_hit.to_bed_fmt(),
        }
    }

    fn to_filter_fmt(&self) -> Vec<String> {
        match self {
            FormattableLine::Annotation(annotation) => annotation.to_filter_fmt(),
            FormattableLine::BenchMarkAnnotation(benchmark) => benchmark.to_filter_fmt(),
            FormattableLine::MaskHit(mask_hit) => mask_hit.to_filter_fmt(),
        }
    }

    fn get_acc(&self) -> String {
        match self {
            FormattableLine::Annotation(annotation) => annotation.get_acc(),
            FormattableLine::BenchMarkAnnotation(benchmark) => benchmark.get_acc(),
            FormattableLine::MaskHit(mask_hit) => mask_hit.get_acc(),
        }
    }
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
        panic!("Data \"{}\" Does Not Exist", assembly_path);
    }
    let fam_file: String = format!("{}/{}.bed.bgz", &assembly_path, &fam);
    if !Path::new(&fam_file).exists() {
        panic!("Family {} Not Found In Assembly {}", &fam, assembly_path);
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
        header = "#seq_id\tseq_start\tseq_end\tfamily_accession\thit_bit_score\tstrand\tbias\tali_start\tali_end\tmodel_start\tmodel_end\thit_evalue_score\tnrph_hit\tdivergence\t*family_name\tseq_len\t*cigar\t*caf";
    }

    writer
        .write_all(format!("{}\n", header).as_bytes())
        .expect("Unable to write line");

    let mut hmm_len = "0".to_string();
    let mut seq_data: Value = json!(null);
    if dl_fmt {
        hmm_len = match json_query(
            &assembly,
            &MOD_LEN_DIR.to_string(),
            &fam.to_string(),
            &"length".to_string(),
            Some(data_dir),
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
        seq_data = seq_json
            .get("data")
            .expect(&format!("Sequence Info For {} Not Found", &fam))
            .to_owned();
    }

    let mut seq_name;
    let mut output;
    for result in reader.lines() {
        let line = result?;
        let fields: Vec<_> = line.split_whitespace().collect();
        let formatted_line = FormattableLine::from_bed(&fields, data_type);
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
                output = formatted_line.to_dl_fmt(seq_name, &hmm_len);
            } else {
                // output = fields.drain(..16).map(|f| f.to_string()).collect(); // TODO readjust this!
                output = formatted_line.to_filter_fmt()
            }
            writer
                .write_all(format!("{}\n", &output.join("\t")).as_bytes())
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
        panic!("Input TSV \"{}\" Not Found", &in_tsv);
    }

    let db_dir = format!("{}/{}", &data_dir, &assembly);
    let target_dir = format!("{}/{}", &db_dir, &data_type);
    if !Path::new(&db_dir).exists() {
        create_dir_all(&db_dir)?;
    }

    if !Path::new(&target_dir).exists() {
        create_dir_all(&target_dir)?;
    }

    let worker_count: NonZeroUsize = match NonZeroUsize::new(10) {
        Some(n) => n,
        None => unreachable!(),
    };

    let in_f = File::open(in_tsv).expect("Could Not Open Input File");
    let lines = BufReader::new(in_f).lines();
    let mut current_acc = "".to_string();
    let mut out_f = tempfile()?;
    let mut out_writer = bgzf::MultithreadedWriter::with_worker_count(worker_count, out_f);
    let mut seen_accs = Vec::new();
    for result in lines {
        let line = result?;
        if !&line.starts_with('#') {
            let fields: Vec<_> = line.split_whitespace().collect();
            let output = FormattableLine::from_export_tsv(&fields, data_type);
            let out_acc = output.get_acc();
            if out_acc != current_acc {
                if seen_accs.contains(&out_acc) {
                    panic!("Input TSV is not sorted in accession order")
                } else {
                    seen_accs.push(out_acc.clone());
                }
                // assume accession order TODO confirm this
                println!("\t{out_acc}");
                current_acc = out_acc;
                out_f = File::create(format!("{target_dir}/{current_acc}.bed.bgz",))
                    .expect("Could Not Open Output File");
                out_writer = bgzf::MultithreadedWriter::with_worker_count(worker_count, out_f);
            };
            out_writer
                .write_all(format!("{}\n", output.to_bed_fmt().join("\t")).as_bytes())
                .expect("Unable to write line");
        }
    }

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
        println!("\tQueued {}: {}", element, needed);
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
                copy(
                    source,
                    format!(
                        "{}/{}{}",
                        target,
                        assembly,
                        file_to_source(element).unwrap()
                    ),
                )
                .expect("Could Not Copy JSON");
                println!("   {} Prep Complete", element);
            } else if source.ends_with(".tsv") {
                println!("   Splitting And Compressing BED Files For {}", element);
                prep_beds(assembly, source, &element.to_string(), Some(data_dir))
                    .expect("BED File Prep Failed");
                println!("   Indexing {}", element);
                let (filenames, bgz_dir, mut contig_index, index_file) =
                    idx::prep_idx(&working_dir, &element.to_string()).expect("Index Prep Failed");
                idx::build_idx(&filenames, &bgz_dir, &mut contig_index, &index_file)
                    .expect("Indexing Failed");
                println!("   {} Prep Complete", element);
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
) -> Result<()> {
    let data_dir = data_directory.unwrap_or(DATA_DIR);
    let assembly_path: String = format!("{}/{}", &data_dir, &assembly_id);
    if !Path::new(&assembly_path).exists() {
        panic!("Assembly \"{}\" Does Not Exist", assembly_path);
    }
    let fam_file: String = format!("{}/{}/{}.bed.bgz", &assembly_path, &ASSEMBLY_DIR, &id);
    if !Path::new(&fam_file).exists() {
        panic!("Family {} Not Found In Assembly {}", id, assembly_path);
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
        Some(data_dir),
    ) {
        Ok(()) => return Ok(()),
        Err(err) => {
            panic!("Error Filtering File: {} - {}", fam_file, err);
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
) -> Result<String> {
    let data_dir = data_directory.unwrap_or(DATA_DIR);
    let assembly_path: String = format!("{}/{}", &data_dir, &assembly);
    // confirm assembly_id and ensure that it accessable
    if !Path::new(&assembly_path).exists() {
        panic!("Assembly \"{}\" Does Not Exist", assembly_path);
    }

    let (filenames, bgz_dir, mut contig_index, index_file) =
        match idx::prep_idx(&assembly_path, &data_type) {
            Ok(res) => res,
            Err(e) => panic!("Search Prep Failed, Index may not exist - {:?}", e),
        };

    if !Path::new(&index_file).exists() {
        panic!(
            "Assembly \"{}\" Is Not Indexed For {}",
            assembly_path, &data_type
        );
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

    let mut formatted = Vec::new();
    match &results {
        Err(e) => {
            panic!("Index Search Failed - {}", e);
        }
        Ok(l) if l.as_slice().is_empty() => {
            return Ok("[]".to_string());
        }
        Ok(l) => {
            for i in 0..l.len() {
                let fields = l[i].split_whitespace().collect::<Vec<&str>>();
                formatted.push(FormattableLine::from_bed(&fields, data_type).to_json());
            }
        }
    };
    match serde_json::to_string(&formatted) {
        Err(e) => {
            panic!("Error Converting Results to JSON - {e}");
        }
        Ok(json_str) => {
            return Ok(json_str);
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
        "{}/{}/{}/{}-{}.json",
        &data_dir, &assembly, &data_type, &assembly, &data_type
    );
    if !Path::new(&target_file).exists() {
        panic!("{} Not Found", &target_file);
    }

    let in_str = read_to_string(&target_file).expect("Could Not Read String");
    let in_data: Value = serde_json::from_str(&in_str).expect("JSON was not well-formatted");

    let val = in_data
        .get("data")
        .and_then(|data| data.get(key).and_then(|item| item.get(target)))
        .expect("Key Target Pair Not Found");
    return Ok(val.to_string().replace("\"", ""));
}

pub fn get_chrom_id(assembly: &String, query: &String) -> String { // TODO write test
    let target_file = format!(
        "{}/{}/{}/{}-{}.json",
        &DATA_DIR, &assembly, "sequences",  &assembly, "sequences"
    );
    if !Path::new(&target_file).exists() {
        eprintln!("{} Not Found", &target_file);
        std::process::exit(1)
    }
    let in_str = read_to_string(&target_file).expect("Could Not Read String");
    let in_json: Value = serde_json::from_str(&in_str).expect("JSON was not well-formatted");
    let data = in_json.get("data").expect("No Data");
    if let Some(data) = data.as_object(){
        for (acc, val) in data {
            if let Some(vals) = val.as_object() {
                let id = vals.get("id").expect("Oh no"); 
                if &id.to_string().replace("\"", "") == query {
                    return acc.to_string();
                }
            }
        }
        panic!("Sequence ID Not Found")
    } else {
        panic!("Something Didn't Work")
    }
}

// OLD Methods ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// pub fn process_json(in_file: &String, key: &String, outfile: &Option<String>) -> Result<()> {
//     if !Path::new(&in_file).exists() {
//         eprintln!("{} Not Found", &in_file);
//         std::process::exit(1)
//     }
//     let in_f = File::open(in_file).expect("Could Not Open Input File");

//     let in_data: Value = serde_json::from_reader(in_f)?;

//     let mut out_data = HashMap::new();
//     if let Some(items) = in_data[2]["data"].as_array() {
//         for item in items {
//             let mut new_item = item.clone().as_object().unwrap().to_owned();
//             new_item.remove(key);
//             out_data.insert(item[key].to_string().replace("\"", ""), new_item);
//         }
//     }
//     let json = json!({
//         "assembly": in_data[1]["name"],
//         "version": in_data[0]["version"],
//         "data": out_data
//     });

//     let output = serde_json::to_string(&json)?;
//     let mut writer: Box<dyn Write> = match outfile {
//         Some(outfile) => Box::new(File::create(&outfile).expect("Could Not Open Output File")),
//         None => Box::new(stdout()),
//     };
//     writer
//         .write_all(format!("{}\n", output).as_bytes())
//         .expect("Unable to write JSON");
//     Ok(())
// }
