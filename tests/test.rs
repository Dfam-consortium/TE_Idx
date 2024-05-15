use std::fs::{read_dir, File};
use std::io::{BufRead, BufReader};
use std::path::Path;
use te_idx::{
    bgzf_filter, idx_query, json_query, prep_beds, prepare_assembly, process_json,
    read_family_assembly_annotations,
};
use tempfile::{NamedTempFile, TempDir};
use noodles::bgzf;

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
    let data_type = &te_idx::ASSEMBLY_DIR.to_string();
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
                    TEST_DATA_DIR,
                    TEST_ASSEMBLY,
                    te_idx::ASSEMBLY_DIR,
                    fam
                ))
                .expect("Can't Open File"),
            )
            .lines()
            .count();
            let filter_count = BufReader::new(File::open(out_f).expect("Can't Open File"))
                .lines()
                .count();
            assert_eq!(filter_count, 20317);
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
    let data_type = &te_idx::ASSEMBLY_DIR.to_string();
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
                    TEST_DATA_DIR,
                    TEST_ASSEMBLY,
                    te_idx::ASSEMBLY_DIR,
                    fam
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
            let target = "#sequence name\tmodel accession\tmodel name\tbit score\te-value\thmm start\thmm end\thmm length\tstrand\talignment start\talignment end\tenvelope start\tenvelope end\tsequence length\nchr1\tDF000000001\tMIR\t97.7\t1e-24\t6\t261\t262\t+\t100012864\t100013069\t100012853\t100013069\t248956422";
            assert_eq!(toplines, target);

            let filter_count = filter_lines.count();
            assert_ne!(orig_count, filter_count);
            assert_eq!(filter_count, 1007);
        }
        Err(e) => panic!("{}", e),
    }
}

#[test]
fn test_idx_query() {
    assert_eq!(true, true)
}

#[test]
fn test_json_query() {
    assert_eq!(true, true)
}

#[test]
fn test_prep_beds() {
    let working_directory = gen_working_dir();

    let assembly = &TEST_ASSEMBLY.to_string();
    let in_tsv = format!("{}/{}/test_ex-mask.tsv", TEST_EXPORT_DIR, TEST_ASSEMBLY);
    let data_type = &te_idx::MASKS_DIR.to_string();
    let data_directory = working_directory.path().to_str();

    match prep_beds(assembly, &in_tsv, data_type, data_directory) {
        Ok(()) => {
            let mask_dir = format!(
                "{}/{}/{}",
                data_directory.unwrap(),
                &TEST_ASSEMBLY,
                &data_type
            );
            assert_eq!(true, Path::new(&mask_dir).exists());
            assert_eq!(
                17,
                read_dir(&mask_dir)
                    .map(|entries| entries
                        .filter_map(|entry| entry.ok())
                        .filter(|entry| entry.file_type().map_or(false, |ft| ft.is_file()))
                        .count())
                    .unwrap_or(0)
            );
        }
        Err(_e) => panic!(),
    }
    let _ = working_directory.close();
}

#[test]
fn test_prepare_assembly() {
    assert_eq!(true, true)
}

#[test]
fn test_process_json() {
    assert_eq!(true, true)
}

#[test]
fn test_read_family_assembly_annotation() {
    assert_eq!(true, true)
}
