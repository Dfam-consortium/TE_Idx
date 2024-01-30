use noodles::bgzf;
use std::fs::{read_dir, File};
use std::io::{stdout, BufRead, Result, Write};
use std::num::NonZeroUsize;
use std::path::Path;
use std::process::exit;

const DATA_DIR: &str = "./data";
const ASSEMBLY_DIR: &str = "assembly_alignments";
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

pub fn read_family_assembly_annotations(id: &u64, assembly_id: &String, nrph: &bool) {
    let assembly_path: String = format!("{}/{}", &DATA_DIR, &assembly_id);
    if !Path::new(&assembly_path).exists() {
        eprintln!("Assembly \"{}\" Does Not Exist", assembly_path);
        exit(1)
    }
    let fam_file: String = format!("{}/{}.bed.bgz", &assembly_path, &id);
    if !Path::new(&fam_file).exists() {
        eprintln!("Family {} Not Found In Assembly {}", id, assembly_path);
        exit(1)
    }

    // if *nrph {
    //     let outfile = format!("{}/{}.nrph.bed.bgz", &assembly_path, &id);
    //     let position: usize = 13;
    //     match bgzf_filter(&fam_file, &position,&String::from("1"), &outfile){
    //         Ok(()) => exit(0),
    //         Err(err) => {eprintln!("Error Filtering File: {} - {}", fam_file, err); exit(1)}
    //     }
    // } else {
    //     exit(0)
    // }
}

pub fn read_annotations(
    assembly: &String,
    chrom: &u32,
    start: &u64,
    end: &u64,
    family: &Option<String>,
    nrph: &bool,
) {
    let assembly_path: String = format!("{}/{}", &DATA_DIR, &assembly);
    // confirm assembly_id and ensure that it accessable
    if !Path::new(&assembly_path).exists() {
        eprintln!("Assembly \"{}\" Does Not Exist", assembly_path);
        exit(1)
    }

    println!(
        "{}, {}, {}, {} ,{}",
        chrom,
        start,
        end,
        if family.is_some() {
            family.as_ref().unwrap()
        } else {
            ""
        },
        nrph
    )
    // find all annotations where either the start or end point is in between the 'start' and 'end' of the window.
    // attributes: ["family_accession", "seq_start", "seq_end", "strand", "ali_start", "ali_end", "model_start", "model_end", "hit_bit_score", "hit_evalue_score", "nrph_hit" sequenceModel.id],
    // if family_accession push where '$family_accession$': family_accession
    // if nrph push nrph_hit: 1

    // collect nhmmerResults Accumulate accessions whose names and types we need to retrieve
    // Retrieve the names and types of all matched families
    // query trfs

    // resolve response {offset: start, length: Math.abs(end - start), query: `${chrom}:${start}-${end}`, hits: nhmmerResults, tandem_repeats: trfResults} 200
}

pub fn find_family(_id: &String, assembly: &String) {
    let paths = read_dir(format!("{}/{}/{}", DATA_DIR, assembly, ASSEMBLY_DIR)).unwrap();
    for path in paths {
        // TODO
        // if path.ends_with(format!("{}.bed.bgz", id)) {
        println!("Name: {}", path.unwrap().path().display())
        // }
    }
}
