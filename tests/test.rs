use noodles::bgzf;
use serde_json::from_str;
use std::collections::HashMap;
use std::fs::{read_dir, File};
use std::io::{BufRead, BufReader};
use std::path::Path;
use te_idx::idx::{build_idx, prep_idx};
use te_idx::{
    bgzf_filter, idx_query, json_query, prep_beds, prepare_assembly,
    read_family_assembly_annotations, ASSEMBLY_DIR, BENCHMARK_DIR, MASKS_DIR, MOD_LEN_DIR,
    SEQUENCE_DIR,
};
use tempfile::{NamedTempFile, TempDir};

pub const TEST_DIR: &'static str = "/home/agray/te_idx/tests";
pub const TEST_DATA_DIR: &'static str = "/home/agray/te_idx/tests/data";
pub const TEST_EXPORT_DIR: &'static str = "/home/agray/te_idx/tests/test_exports";

pub const TEST_ASSEMBLY: &'static str = "test_ex";

fn gen_working_dir() -> TempDir {
    return TempDir::new_in(TEST_DIR).expect("Error Creating Working Directory");
}

#[test]
fn test_bgzf_filter_nrph() {
    let out_f = NamedTempFile::new_in(TEST_DATA_DIR).expect("Couldn't Open Output File");
    let f_name = out_f.path().to_str().unwrap();

    let assembly = &TEST_ASSEMBLY.to_string();
    let data_type = &ASSEMBLY_DIR.to_string();
    let fam = &"DF000000001".to_string();
    let position: usize = 13;
    let term: Option<String> = Some("1".to_string());
    let outfile = Some(f_name.to_string());
    let dl_fmt = false;
    let data_directory = Some(TEST_DATA_DIR);

    // Test for NRPH filter
    match bgzf_filter(
        assembly,
        data_type,
        fam,
        &position,
        &term,
        &outfile,
        dl_fmt,
        data_directory,
    ) {
        Ok(()) => {
            let orig_count = BufReader::new(
                File::open(format!(
                    "{}/{}/{}/{}.bed.bgz",
                    TEST_DATA_DIR, TEST_ASSEMBLY, ASSEMBLY_DIR, fam
                ))
                .expect("Can't Open File"),
            )
            .lines()
            .count();
            let filter_count = BufReader::new(File::open(out_f).expect("Can't Open File"))
                .lines()
                .count();
            // check that filtered file is smaller and contains expected lines
            assert_eq!(filter_count, 5975);
            assert_ne!(orig_count, filter_count);
        }
        Err(e) => panic!("{}", e),
    }
}

#[test]
fn test_bgzf_filter_dl_fmt() {
    let out_f = NamedTempFile::new_in(TEST_DATA_DIR).expect("Couldn't Open Output File");
    let f_name = out_f.path().to_str().unwrap();

    let assembly = &TEST_ASSEMBLY.to_string();
    let data_type = &ASSEMBLY_DIR.to_string();
    let fam = &"DF000000001".to_string();
    let position = 7;
    let term = Some("14.7".to_string());
    let outfile = Some(f_name.to_string());
    let dl_fmt = true;
    let data_directory = Some(TEST_DATA_DIR);

    // Test for download format
    match bgzf_filter(
        assembly,
        data_type,
        fam,
        &position,
        &term,
        &outfile,
        dl_fmt,
        data_directory,
    ) {
        Ok(()) => {
            let orig_count = BufReader::new(
                File::open(format!(
                    "{}/{}/{}/{}.bed.bgz",
                    TEST_DATA_DIR, TEST_ASSEMBLY, ASSEMBLY_DIR, fam
                ))
                .expect("Can't Open File"),
            )
            .lines()
            .count();
            let mut filter_lines =
                bgzf::Reader::new(File::open(out_f).expect("Can't Open File")).lines();

            let first = filter_lines
                .next()
                .unwrap_or_else(|| Ok("".to_string()))
                .expect("No Line");
            let second = filter_lines
                .next()
                .unwrap_or_else(|| Ok("".to_string()))
                .expect("No Line");
            let toplines = format!("{}\n{}", first, second);
            let target = "#sequence name\tmodel accession\tmodel name\tbit score\te-value\thmm start\thmm end\thmm length\tstrand\talignment start\talignment end\tenvelope start\tenvelope end\tsequence length\nchr1\tDF000000001\tMIR\t11.0\t290.0\t37\t206\t262\t-\t89176700\t89176576\t89176720\t89176556\t248956422";
            // check that header and first line are correct
            assert_eq!(toplines, target);

            // check that filtered file is smaller and contains expected lines
            let filter_count = filter_lines.count();
            assert_ne!(orig_count, filter_count);
            assert_eq!(filter_count, 96);
        }
        Err(e) => panic!("{}", e),
    }
}

#[test]
fn test_build_idx() {
    let data_dir = TEST_DATA_DIR;
    let assembly = TEST_ASSEMBLY;
    let data_type = &MASKS_DIR.to_string();

    let (filenames, bgz_dir, mut contig_index, index_file) =
        match prep_idx(&format!("{}/{}", &data_dir, &assembly), data_type) {
            Ok(res) => res,
            Err(e) => panic!(
                "Search Prep Failed, Assembly or Data Type May Not Exist - {:?}",
                e
            ),
        };
    assert_eq!(filenames.len(), 341);
    assert_eq!(bgz_dir, "/home/agray/te_idx/tests/data/test_ex/masks");
    assert_eq!(
        index_file,
        "/home/agray/te_idx/tests/data/test_ex/masks_idx.dat"
    );

    let working_directory = gen_working_dir();
    let test_data_dir = working_directory.path().to_str().unwrap();
    let test_index_file = format!("{}/masks_idx.dat", test_data_dir);

    build_idx(&filenames, &bgz_dir, &mut contig_index, &test_index_file).expect("Indexing Failed");

    assert!(Path::new(&test_index_file).exists());

    let _ = working_directory.close();
}

#[test]
fn test_idx_query() {
    let assembly = &TEST_ASSEMBLY.to_string();
    let data_type = &ASSEMBLY_DIR.to_string();
    let chrom = &1.to_string();
    let start = 10000;
    let end = 100000;
    let family: &Option<String> = &None;
    let nrph = &false;
    let data_directory = Some(TEST_DATA_DIR);

    let res1 = idx_query(
        assembly,
        data_type,
        chrom,
        start,
        end,
        family,
        nrph,
        data_directory,
    )
    .expect("Index Query Failed");
    let vals1: Vec<HashMap<String, String>> = from_str(&res1).expect("Cannot Deserialize");
    assert_eq!(vals1.len(), 22);
}

#[test]
fn test_idx_query_fam() {
    let assembly = &TEST_ASSEMBLY.to_string();
    let data_type = &ASSEMBLY_DIR.to_string();
    let chrom = &1.to_string();
    let start = 10000;
    let end = 100000;
    let family: &Option<String> = &Some("DF000000001".to_string());
    let nrph = &false;
    let data_directory = Some(TEST_DATA_DIR);

    let res = idx_query(
        assembly,
        data_type,
        chrom,
        start,
        end,
        family,
        nrph,
        data_directory,
    )
    .expect("Index Query Failed");
    println!("{:?}", res);
    let vals: Vec<HashMap<String, String>> = from_str(&res).expect("Cannot Deserialize");
    assert_eq!(vals.len(), 15);
}

#[test]
fn test_idx_query_nrph() {
    let assembly = &TEST_ASSEMBLY.to_string();
    let data_type = &ASSEMBLY_DIR.to_string();
    let chrom = &1.to_string();
    let start = 10000;
    let end = 100000;
    let family: &Option<String> = &Some("DF000000001".to_string());

    let nrph = &true;
    let data_directory = Some(TEST_DATA_DIR);

    let res = idx_query(
        assembly,
        data_type,
        chrom,
        start,
        end,
        family,
        nrph,
        data_directory,
    )
    .expect("Index Query Failed");
    let vals: Vec<HashMap<String, String>> = from_str(&res).expect("Cannot Deserialize");
    assert_eq!(vals.len(), 15);
}

#[test]
fn test_json_query() {
    let assembly = &TEST_ASSEMBLY.to_string();
    let data_type = &SEQUENCE_DIR.to_string();
    let key = &"1".to_string();
    let target = &"length".to_string();
    let data_directory = Some(TEST_DATA_DIR);

    let ans =
        json_query(assembly, data_type, key, target, data_directory).expect("JSON Read Failed");
    assert_eq!(ans, "248956422");
}

#[test]
fn test_prep_beds() {
    let working_directory = gen_working_dir();

    let assembly = &TEST_ASSEMBLY.to_string();
    let in_tsv = format!("{}/{}/test_ex-mask.tsv", TEST_EXPORT_DIR, TEST_ASSEMBLY);
    let data_type = &MASKS_DIR.to_string();
    let data_directory = working_directory.path().to_str();

    match prep_beds(assembly, &in_tsv, data_type, data_directory) {
        Ok(()) => {
            let mask_dir = format!(
                "{}/{}/{}",
                data_directory.unwrap(),
                &TEST_ASSEMBLY,
                &data_type
            );
            // check that new folder was created and contains expected number of files
            assert_eq!(true, Path::new(&mask_dir).exists());
            assert_eq!(
                19,
                read_dir(&mask_dir)
                    .map(|entries| entries
                        .filter_map(|entry| entry.ok())
                        .filter(|entry| entry.file_type().map_or(false, |ft| ft.is_file()))
                        .count())
                    .unwrap_or(0)
            );
        }
        Err(_e) => {
            let _ = working_directory.close();
            panic!()
        }
    }
    let _ = working_directory.close();
}

#[test]
fn test_prepare_assembly() {
    let working_directory = gen_working_dir();
    let test_data_dir = Some(working_directory.path().to_str().unwrap());

    let assembly = &TEST_ASSEMBLY.to_string();
    let data_directory = test_data_dir;
    let export_directory = Some(TEST_EXPORT_DIR);

    let _ = prepare_assembly(assembly, data_directory, export_directory);
    let assembly_dir = format!("{}/{}", data_directory.unwrap(), TEST_ASSEMBLY);

    let align_dir = &format!("{}/{}", assembly_dir, ASSEMBLY_DIR);
    let mut align_contents = Path::new(&align_dir).read_dir().expect("Couldn't Read Dir");
    let bench_dir = &format!("{}/{}", assembly_dir, BENCHMARK_DIR);
    let mut bench_contents = Path::new(&bench_dir).read_dir().expect("Couldn't Read Dir");
    let mask_dir = &format!("{}/{}", assembly_dir, MASKS_DIR);
    let mut mask_contents = Path::new(&mask_dir).read_dir().expect("Couldn't Read Dir");
    let mlen_dir = &format!("{}/{}", assembly_dir, MOD_LEN_DIR);
    let mut mlen_contents = Path::new(&mlen_dir).read_dir().expect("Couldn't Read Dir");
    let seq_dir = &format!("{}/{}", assembly_dir, SEQUENCE_DIR);
    let mut seq_contents = Path::new(&seq_dir).read_dir().expect("Couldn't Read Dir");

    assert!(Path::new(align_dir).exists());
    assert!(Path::new(bench_dir).exists());
    assert!(Path::new(mask_dir).exists());
    assert!(Path::new(mlen_dir).exists());
    assert!(Path::new(seq_dir).exists());
    assert!(Path::new(&format!("{}/{}_idx.dat", assembly_dir, ASSEMBLY_DIR)).exists());
    assert!(Path::new(&format!("{}/{}_idx.dat", assembly_dir, BENCHMARK_DIR)).exists());
    assert!(Path::new(&format!("{}/{}_idx.dat", assembly_dir, MASKS_DIR)).exists());
    assert!(!align_contents.next().is_none());
    assert!(!bench_contents.next().is_none());
    assert!(!mask_contents.next().is_none());
    assert!(!mlen_contents.next().is_none());
    assert!(!seq_contents.next().is_none());

    let _c = working_directory.close();
}

#[test]
fn test_read_family_assembly_annotation() {
    let out_f = NamedTempFile::new_in(TEST_DATA_DIR).expect("Couldn't Open Output File");
    let id = &"DF000000001".to_string();
    let assembly_id = &TEST_ASSEMBLY.to_string();
    let nrph = &false;
    let outfile = &Some(out_f.path().to_str().unwrap().to_string());
    let data_directory = Some(TEST_DATA_DIR);

    let _ = read_family_assembly_annotations(id, assembly_id, nrph, outfile, data_directory);

    let reader = bgzf::Reader::new(File::open(out_f).expect("can't open"));
    let line_count = reader.lines().count();
    assert_eq!(line_count, 60001);
}

#[test]
fn test_read_family_assembly_annotation_nrph() {
    let out_f = NamedTempFile::new_in(TEST_DATA_DIR).expect("Couldn't Open Output File");

    let id = &"DF000000001".to_string();
    let assembly_id = &TEST_ASSEMBLY.to_string();
    let nrph = &true;
    let outfile = &Some(out_f.path().to_str().unwrap().to_string());
    let data_directory = Some(TEST_DATA_DIR);

    let _ = read_family_assembly_annotations(id, assembly_id, nrph, outfile, data_directory);

    let reader = bgzf::Reader::new(File::open(out_f).expect("can't open"));
    let line_count = reader.lines().count();
    assert_eq!(line_count, 59911);
}

// OLD
// #[test]
// fn test_process_json() {
//     let out_f = NamedTempFile::new_in(TEST_DATA_DIR).expect("Couldn't Open Output File");

//     let in_file = &format!(
//         "{}/{}/{}-{}.json",
//         TEST_EXPORT_DIR, TEST_ASSEMBLY, TEST_ASSEMBLY, SEQUENCE_DIR
//     );
//     let key = &"accession".to_string();
//     let outfile = &Some(out_f.path().to_str().unwrap().to_string());

//     process_json(in_file, key, outfile).expect("JSON Parse Failed");

//     let output = read_to_string(out_f).expect("Can't Read File");
//     let comparison = format!(
//         "{}/{}/{}/{}-processed-{}.json",
//         TEST_DATA_DIR, TEST_ASSEMBLY, SEQUENCE_DIR, TEST_ASSEMBLY, SEQUENCE_DIR
//     );
//     let target = read_to_string(comparison).expect("Can't Read File");
//     assert_eq!(output, target);
// }
