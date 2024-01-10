use clap::Parser;
use std::io::{Result, Write, BufRead};
use std::path::Path;
use std::fs::File;
use noodles::bgzf as bgzf;
use std::num::NonZeroUsize;

#[derive(Parser)]
#[derive(Debug)]
struct Cli {
    #[arg(long, short)]
    infile: String,
    #[arg(long, short)]
    position: usize,
    #[arg(long, short)]
    term: String,
    #[arg(long, short)]
    outfile: String
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    if !Path::new(&cli.infile).exists() {
        println!("{} Not Found", cli.infile);
        std::process::exit(1)
    }

    let worker_count: NonZeroUsize = match NonZeroUsize::new(5) {
        Some(n) => n,
        None => unreachable!(),
    };
    let in_f = File::open(cli.infile).expect("Could Not Open Input File");
    let out_f = File::create(cli.outfile).expect("Could Not Open Output File");

    let reader = bgzf::MultithreadedReader::with_worker_count(worker_count, in_f);
    let mut writer = bgzf::MultithreadedWriter::with_worker_count(worker_count, out_f);

    for line in reader.lines() {
        let linestr = line.unwrap();
        let mut fields: Vec<_> = linestr.split_whitespace().collect();
        if fields.get(cli.position - 1).unwrap() == &&cli.term {
            fields.truncate(9);
            writer.write_all(format!("{}\n", fields.join("\t")).as_bytes()).expect("Unable to write line");
        }
    }
    Ok(())
}
