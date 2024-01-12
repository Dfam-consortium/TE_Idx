use noodles::bgzf;
use std::fs::{read_dir, File};
use std::io::{BufRead, Result, Write};
use std::num::NonZeroUsize;
use std::path::Path;

const DATA_DIR: &str = "./data";
const ASSEMBLY_DIR: &str = "assembly_alignments";
// const BENCHMARK_DIR: &str = "./data/benchmark_alignments";
// const SEQUENCE_DIR: &str = "./data/sequence";

pub fn bgzf_filter(
    infile: &String,
    position: &usize,
    term: &String,
    outfile: &String,
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
    let out_f = File::create(outfile).expect("Could Not Open Output File");

    let reader = bgzf::MultithreadedReader::with_worker_count(worker_count, in_f);
    let mut writer = bgzf::MultithreadedWriter::with_worker_count(worker_count, out_f);

    for line in reader.lines() {
        let linestr = line.unwrap();
        let mut fields: Vec<_> = linestr.split_whitespace().collect();
        if fields.get(position - 1).unwrap() == term {
            fields.truncate(9);
            writer
                .write_all(format!("{}\n", fields.join("\t")).as_bytes())
                .expect("Unable to write line");
        }
    }
    Ok(())
}

pub fn read_family_assembly_annotations(
    _id: &String,
    _assembly_id: &String,
    nrph: bool,
    _download: bool,
) {
    // Retrieve all annotations

    // id -> family id (DF000000001)
    // assembly_id -> assembly (hg38)
    // nrph -> bool
    // download -> bool

    // confirm assembly_id and ensure that it accessable
    // if not, return 404
    let _column = if nrph { "nrph_hit_list" } else { "hit_list" };
    // query model_file for column with id
    // if nothing found, return 404
    // return 200, gzip text of data values
}

pub fn read_annotations(
    _assembly: &String,
    _chrom: &usize,
    _start: &usize,
    _end: &usize,
    _family: &String,
    _nrph: bool,
) {
    // retrieve all annotations from range in assembly

    // family_accession = family;
    // if (start > end) {
    //     let swap = start;
    //     start = end;
    //     end = swap;
    // }

    // |end-start| > 1000000 return 404  message: "Requested range is too long."

    // confirm assembly_id and ensure that it accessable
    // if not, return 404

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

        // if path.ends_with(format!("{}.bed.bgz", id)) {
            println!("Name: {}", path.unwrap().path().display())
        // }
    }
}
