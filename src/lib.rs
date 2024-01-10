use std::io::{Result, Write, BufRead};
use std::path::Path;
use std::fs::{File, read_dir};
use noodles::bgzf as bgzf;
use std::num::NonZeroUsize;

const ASSEMBLY_DIR: &str = "./data/assembly_alignments";
// const BENCHMARK_DIR: &str = "./data/benchmark_alignments";
// const SEQUENCE_DIR: &str = "./data/sequence";

pub fn bgzf_filter(infile: &String, position: &usize, term: &String, outfile: &String)  -> Result<()>{
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
            writer.write_all(format!("{}\n", fields.join("\t")).as_bytes()).expect("Unable to write line");
        }
    }
    Ok(())
}

pub fn read_family_assembly_annotations(_id: String, _assembly_id: String, nrph: bool, _download: bool){
    // id -> family id (DF000000001)
    // assembly_id -> assembly (hg38)
    // nrph -> bool
    // download -> bool

    // confirm assembly_id and ensure that it accessable
    // if not, return 404
    let _column = if nrph {"nrph_hit_list"} else {"hit_list"};
    // query model_file for column with id
    // if nothing found, return 404
    // return 200, gzip text of data values
}

pub fn find_family(_id: &String){
    let paths = read_dir(ASSEMBLY_DIR).unwrap();
    for path in paths {
        println!("Name: {}", path.unwrap().path().display())
    }
}