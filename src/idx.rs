use byteorder::{LittleEndian, ReadBytesExt};
use log::{debug, error, info, warn, Level, LevelFilter, Metadata, Record};
use noodles::bgzf;
use std::collections::HashMap;
use std::collections::HashSet;
use std::error::Error;
use std::fs;
use std::fs::File;
use std::io;
use std::io::prelude::*;
use std::io::{Seek, SeekFrom};
use std::time::SystemTime;

static MY_LOGGER: MyLogger = MyLogger;

struct MyLogger;

impl log::Log for MyLogger {
    fn enabled(&self, metadata: &Metadata) -> bool {
        metadata.level() <= Level::Info
    }

    fn log(&self, record: &Record) {
        if self.enabled(record.metadata()) {
            println!("{} - {}", record.level(), record.args());
        }
    }
    fn flush(&self) {}
}
// Magic number for this index format (6-bytes)
const MAGIC_NUMBER: &[u8] = b"#R_IDX";

// File format version (2-bytes)
const FORMAT_VERSION: u16 = 0;

#[derive(Debug)]
pub struct ContigIndex {
    tile_size: u32, // Default: 16384
    bgz_files: Vec<BGZFile>,
    contig_count: u32,
    tile_counts: Vec<u32>,
    range_counts: Vec<Vec<u32>>,
    range_data_index: Vec<Vec<u32>>,
    contigs: Vec<Contig>,
    contig_lookup: HashMap<String, u32>,
}

#[derive(Debug)]
struct BGZFile {
    name: String,
    mod_time: f64,
    bytes: u64,
}

#[derive(Debug)]
struct Contig {
    name: String,
    contig_tiles: Vec<ContigTile>,
}

#[derive(Debug)]
struct ContigTile {
    contig_ranges: Vec<ContigRange>,
}

impl ContigTile {
    fn new() -> ContigTile {
        ContigTile {
            contig_ranges: Vec::new(),
        }
    }
}

#[derive(Debug, Eq, Ord, PartialEq, PartialOrd, Clone)]
#[repr(C)]
pub struct ContigRange {
    bed_idx: u32,  // Identifies the source BED file for this entry
    start_bp: u64, // region start bp (zero-based)
    end_bp: u64,   // region end bp (zero-based, half-open)
    bgzf_pos: u64, // Byte position within bgz compressed BED file
}

fn read_u16_from_file(file: &mut File) -> io::Result<u16> {
    let mut buffer = [0; 2];
    file.read_exact(&mut buffer)?;
    Ok(u16::from_le_bytes(buffer))
}

fn read_u32_from_file(file: &mut File) -> io::Result<u32> {
    let mut buffer = [0; 4];
    file.read_exact(&mut buffer)?;
    Ok(u32::from_le_bytes(buffer))
}

// fn read_u64_from_file(file: &mut File) -> io::Result<u64> {
//     let mut buffer = [0; 8];
//     file.read_exact(&mut buffer)?;
//     Ok(u64::from_le_bytes(buffer))
// }

// fn read_string_from_file(file: &mut File, length: usize) -> io::Result<String> {
//     let mut buffer = vec![0; length];
//     file.read_exact(&mut buffer)?;
//     Ok(String::from_utf8_lossy(&buffer).to_string())
// }

impl ContigIndex {
    fn get_or_insert_contig(&mut self, contig_name: &str) -> &mut Contig {
        // Check if contig name is already defined in the lookup table
        if !self.contig_lookup.contains_key(contig_name) {
            // Missing, add the new contig
            let contig_idx = self.contigs.len() as u32;
            self.contigs.push(Contig {
                name: contig_name.to_string(),
                contig_tiles: Vec::new(),
            });
            self.contig_lookup
                .insert(contig_name.to_string(), contig_idx);
            return &mut self.contigs[contig_idx as usize];
        }
        // Pre-existing, simply return it
        let contig_idx = *self.contig_lookup.get(contig_name).unwrap();

        &mut self.contigs[contig_idx as usize]
    }

    // Add a new range to the ContigIndex
    fn add_contig_range(
        &mut self,
        contig_name: &str,
        bed_idx: u32,
        start_bp: u64,
        end_bp: u64,
        bgzf_pos: u64,
    ) {
        let first_tile_idx = (start_bp / u64::from(self.tile_size)) as usize;
        let last_tile_idx = ((end_bp - 1) / u64::from(self.tile_size)) as usize;
        let contig = self.get_or_insert_contig(contig_name);
        // Initialize any missing tiles up and including last_tile_idx if not already present.
        for _ in 0..((last_tile_idx + 1).saturating_sub(contig.contig_tiles.len())) {
            contig.contig_tiles.push(ContigTile::new());
        }
        for tile_idx in first_tile_idx..(last_tile_idx + 1) {
            if let Some(tile) = contig.contig_tiles.get_mut(tile_idx) {
                let new_contig_range = ContigRange {
                    bed_idx,
                    start_bp,
                    end_bp,
                    bgzf_pos,
                };
                tile.contig_ranges.push(new_contig_range);
            }
        }
    }

    fn init_search(&mut self, file_path: &str) {
        let mut file = match File::open(file_path) {
            Ok(file) => file,
            Err(e) => panic!("Error Opening File - {:?}", e),
        };

        // Read the magic number (6-bytes)
        let mut buffer = [0; 6];
        file.read_exact(&mut buffer).unwrap();
        if buffer != MAGIC_NUMBER {
            warn!("Magic number did not match.");
        }

        // Read the file format version (2-bytes, little-endian)
        let f_ver = read_u16_from_file(&mut file).unwrap();
        if f_ver != FORMAT_VERSION {
            warn!(
                "Incompatible file version {}, expected {}",
                f_ver, FORMAT_VERSION
            );
        }

        // Read the tile_size used in this index (u32, little-endian)
        self.tile_size = read_u32_from_file(&mut file).unwrap();
        info!("Round trip tile size = {}", self.tile_size);

        // Read the contig count in this index (u32, little-endian)
        self.contig_count = read_u32_from_file(&mut file).unwrap();
        info!("Contig count = {}", self.contig_count);

        // Read the file count in this index (u32, little-endian)
        let file_count = read_u32_from_file(&mut file).unwrap();
        info!("File count = {}", file_count);

        // Read the tile counts (tile_counts[contig], u32, little-endian)
        self.tile_counts = vec![0; self.contig_count as usize];
        file.read_u32_into::<LittleEndian>(&mut self.tile_counts)
            .unwrap();

        // The range data starts at: 20+(files*56)+(Contigs*44)+(Tiles*4)
        let mut range_data_start = 20 + (file_count * 56) + (self.contig_count * 44);
        for c in &self.tile_counts {
            range_data_start += c * 4;
        }

        // Read the range counts (range_counts[contig][tile], u32, little-endian)
        let mut tile_start = range_data_start;
        for contig_idx in 0..self.contig_count {
            // build range_counts
            let mut r_counts = vec![0; self.tile_counts[contig_idx as usize] as usize];
            file.read_u32_into::<LittleEndian>(&mut r_counts).unwrap();
            self.range_counts.push(r_counts);

            // build range_data_index[contig][tile] = file_byte_position
            self.range_data_index.push(Vec::new());
            for tile_idx in 0..self.tile_counts[contig_idx as usize] {
                self.range_data_index[contig_idx as usize].push(tile_start);
                tile_start += self.range_counts[contig_idx as usize][tile_idx as usize] * 28;
            }
        }

        // Read the contig name strings into a buffer
        let mut buffer = vec![0; (40 * self.contig_count) as usize];
        file.read_exact(&mut buffer).unwrap();

        // TODO: Also make a hash from this for reverse lookups
        // Process the buffer as 40-character strings
        let mut offset = 0;
        for n in 0..self.contig_count {
            let string_slice = &buffer[offset..offset + 40];
            // Find the index of the null byte
            let null_byte_index = string_slice.iter().position(|&x| x == b'\0');
            // Create a new string up to the null byte index
            let string = match null_byte_index {
                Some(index) => String::from_utf8_lossy(&string_slice[..index]).to_string(),
                None => String::from_utf8_lossy(string_slice).to_string(),
            };

            self.contig_lookup.insert(string.clone(), n as u32);
            offset += 40;
        }

        // Finally read in the filenames and stats
        let mut buffer = vec![0; (56 * file_count) as usize];
        file.read_exact(&mut buffer).unwrap();
        for i in 0..file_count {
            let mut start = (i * 56) as usize;

            let string_slice = &buffer[start..start + 40];
            // Find the index of the null byte
            let null_byte_index = string_slice.iter().position(|&x| x == b'\0');
            // Create a new string up to the null byte index
            let string = match null_byte_index {
                Some(index) => String::from_utf8_lossy(&string_slice[..index]).to_string(),
                None => String::from_utf8_lossy(string_slice).to_string(),
            };

            start += 40;

            let value1 = f64::from_le_bytes([
                buffer[start],
                buffer[start + 1],
                buffer[start + 2],
                buffer[start + 3],
                buffer[start + 4],
                buffer[start + 5],
                buffer[start + 6],
                buffer[start + 7],
            ]);

            let value2 = u64::from_le_bytes([
                buffer[start + 8],
                buffer[start + 9],
                buffer[start + 10],
                buffer[start + 11],
                buffer[start + 12],
                buffer[start + 13],
                buffer[start + 14],
                buffer[start + 15],
            ]);

            let my_struct = BGZFile {
                name: string,
                mod_time: value1,
                bytes: value2,
            };
            self.bgz_files.push(my_struct);
        }
        // NEW: BGZFiles starts at: 20+(Contigs*44)+(Tiles*4)
        // NEW: ContigRanges starts at: 20+(Files*56)+(Contigs*44)+(Tiles*4)
    }

    fn read_tile(
        &mut self,
        file: &mut File,
        contig: u32,
        tile: usize,
    ) -> Result<ContigTile, &'static str> {
        let range_count = self.range_counts[contig as usize][tile as usize];
        let mut buffer = vec![0; (28 * range_count) as usize];
        let byte_pos = self.range_data_index[contig as usize][tile as usize];
        info!(
            "read_tile: expecting {} ranges, byte_pos {}",
            range_count, byte_pos
        );
        let _ = file.seek(SeekFrom::Start(byte_pos as u64));
        file.read_exact(&mut buffer).unwrap();

        let mut c_tile = ContigTile {
            contig_ranges: Vec::with_capacity(range_count as usize),
        };

        for i in 0..range_count {
            let start = (i * 28) as usize;

            let value0 = u32::from_le_bytes([
                buffer[start],
                buffer[start + 1],
                buffer[start + 2],
                buffer[start + 3],
            ]);

            let value1 = u64::from_le_bytes([
                buffer[start + 4],
                buffer[start + 5],
                buffer[start + 6],
                buffer[start + 7],
                buffer[start + 8],
                buffer[start + 9],
                buffer[start + 10],
                buffer[start + 11],
            ]);

            let value2 = u64::from_le_bytes([
                buffer[start + 12],
                buffer[start + 13],
                buffer[start + 14],
                buffer[start + 15],
                buffer[start + 16],
                buffer[start + 17],
                buffer[start + 18],
                buffer[start + 19],
            ]);

            let value3 = u64::from_le_bytes([
                buffer[start + 20],
                buffer[start + 21],
                buffer[start + 22],
                buffer[start + 23],
                buffer[start + 24],
                buffer[start + 25],
                buffer[start + 26],
                buffer[start + 27],
            ]);

            let my_struct = ContigRange {
                bed_idx: value0,
                start_bp: value1,
                end_bp: value2,
                bgzf_pos: value3,
            };
            c_tile.contig_ranges.push(my_struct);
        }

        Ok(c_tile)
    }

    // TODO: deprecate filenames and store in index
    fn search(
        &mut self,
        i_file: &mut File,
        bgz_dir: &String,
        q_contig: &String,
        q_start: u64,
        q_end: u64,
        q_family: &Option<String>,
        q_nrph: bool,
    ) -> Result<Vec<String>, Box<dyn Error>> {
        fn filter_line(line: &String, q_family: &Option<String>, q_nrph: &bool) -> bool {
            let fields = line.split_whitespace().collect::<Vec<&str>>();
            if q_family.is_some() {
                let acc: &str = fields[3].split(".").collect::<Vec<&str>>()[0];
                if q_family.as_ref().unwrap() != acc {
                    return false;
                };
            }
            if *q_nrph == true {
                if fields.last().unwrap() != &"1" {
                    return false;
                }
            }
            return true;
        }

        // TODO: Return if cannot identify contig
        let q_contig_idx = *self.contig_lookup.get(q_contig).unwrap();

        let mut results: Vec<String> = Vec::new();

        // Determine the start/end tiles this range could possibly overlap
        let start_tile = (q_start / self.tile_size as u64) as usize;
        let mut end_tile = ((q_end - 1) / self.tile_size as u64) as usize;

        // The start position is outside the indexed size of
        // this contig.  Being kind and returning early with
        // null results.
        //   TODO: This should probably be an error as with having
        //         an end outside the index.
        if start_tile > self.tile_counts[q_contig_idx as usize] as usize {
            let error_str = "Start Position Outside Of Indexed Size";
            error!("{}", error_str);
            return Err(error_str.into());
        }

        // The end position is outside the indexed size of
        // the contig.  This is also an error however in this
        // case we could possibly return some overlapping
        // ranges since the start position is overlapping.
        end_tile = end_tile.min((self.tile_counts[q_contig_idx as usize] - 1) as usize);

        info!(
            "Query: {}:{}-{}  tiles:{} to {}",
            q_contig, q_start, q_end, start_tile, end_tile
        );
        let mut hits = 0;

        //  This could be adaptive based on how many tiles it needs to access
        //     e.g if its less than 1000 read in a block and if it's greater
        //         than read record by record as they did.

        let mut range_count = self.range_counts[q_contig_idx as usize][start_tile as usize];
        info!("search: range_count {}", range_count);
        if range_count > 0 {
            let range_data = self
                .read_tile(i_file, q_contig_idx, start_tile)
                .unwrap()
                .contig_ranges;
            if range_data[0].start_bp < q_end {
                // TODO: This has to be the right-most variant as we only need to guarantee
                // that the start position is sorted.

                // Binary search (left-most variant) for first range to the right of query (half open query)
                //
                // Searching for the highest start range first and working backwards
                // allows us to avoid having to have the end position sorted as well
                // as the start position.
                let mut left = 0;
                let mut right = range_count;
                while left < right {
                    // Effectively a floor truncation in Rust
                    let mid = (left + right) / 2;
                    // End is a half-open coordinate
                    if range_data[mid as usize].start_bp < q_end {
                        left = mid + 1;
                    } else {
                        right = mid;
                    }
                }
                // Left is not inclusive so we can use directly in Rust range 0..left
                let mut ranges = Vec::new();
                for r_idx in (0..left).rev() {
                    if range_data[r_idx as usize].start_bp < q_end {
                        ranges.push(range_data[r_idx as usize].clone());
                    } else {
                        break;
                    }
                }
                for range in ranges.iter().rev() {
                    // This is surprisingly fast despite having to open/abandon a bgzf file per
                    // annotation.  Pre-grouping the annotations by family/start might speed up
                    // retreival, however then it would need to be resorted by contig/start for
                    // output -- all in memory -- should experiment.
                    let bgz_file = format!(
                        "{}/{}",
                        bgz_dir, self.bgz_files[range.bed_idx as usize].name
                    );
                    let mut reader = File::open(&bgz_file).map(bgzf::Reader::new).unwrap();
                    reader
                        .seek(bgzf::VirtualPosition::from(range.bgzf_pos))
                        .unwrap();
                    let mut line = String::new();
                    reader.read_line(&mut line).unwrap();
                    if filter_line(&line, &q_family, &q_nrph) {
                        results.push(line);
                        hits += 1;
                    }
                }
            }
            if end_tile > start_tile {
                let mut tile_start_bp = (self.tile_size as u64) * ((start_tile + 1) as u64);
                for t_idx in (start_tile + 1)..=end_tile {
                    range_count = self.range_counts[q_contig_idx as usize][t_idx];
                    if range_count > 0 {
                        let range_data = self
                            .read_tile(i_file, q_contig_idx, t_idx)
                            .unwrap()
                            .contig_ranges;
                        if range_data[0].start_bp < q_end {
                            // A binary search is not needed here as we know that the query spans
                            // more than one tile and that either this tile needs to be evaluated
                            // fully (middle tile), or contains hits on the left side and simply
                            // needs a break condition when it has past the last annotation
                            // (end tile).
                            for r_idx in 0..range_count {
                                if range_data[r_idx as usize].start_bp < (tile_start_bp as u64) {
                                    continue;
                                }
                                if range_data[r_idx as usize].start_bp < q_end {
                                    let bgz_file = format!(
                                        "{}/{}",
                                        bgz_dir,
                                        self.bgz_files[range_data[r_idx as usize].bed_idx as usize]
                                            .name
                                    );
                                    let mut reader =
                                        File::open(&bgz_file).map(bgzf::Reader::new).unwrap();
                                    reader
                                        .seek(bgzf::VirtualPosition::from(
                                            range_data[r_idx as usize].bgzf_pos,
                                        ))
                                        .unwrap();
                                    let mut line = String::new();
                                    reader.read_line(&mut line).unwrap();
                                    if filter_line(&line, &q_family, &q_nrph) {
                                        results.push(line);
                                        hits += 1;
                                    }
                                } else {
                                    info!(
                                        "Breaking because {} >= {}",
                                        range_data[r_idx as usize].start_bp, q_end
                                    );
                                    break;
                                }
                            }
                        }
                    }
                    tile_start_bp += self.tile_size as u64;
                }
            }
        }
        info!("Total overlaps: {}", hits);
        Ok(results)
    }

    //
    // Perhaps (niavely) I didn't use Serde to do this. I was
    // worried that it's serialization mechanism was too opaque
    // to support dependable file reading using seek.
    //
    // Index File:
    //  Bytes   Type   Byte_order      Description
    //  -----   -----  --------------  -------------------------
    //   6      chars                  Magic Number "#R_IDX"
    //   2      u16    little-endian   Format Version
    //   4      u32    little-endian   Tile Size
    //   4      u32    little-endian   Contig Count (C)
    //   4      u32    little-endian   BGZ File Count (F) **NEW**
    //  C*4     [u32]  little-endian   Per contig tile counts (T)
    //  T*4     [u32]  little-endian   Per tile range counts (R)
    //  C*40    chars                  Per contig 40 char name
    //  F*56    BGZ_File (see details) **NEW**
    //  R*28    ContigRanges (see details) in contig,tile order
    //
    // BGZ_File Structure (56 bytes)
    //  Bytes   Type   Byte_order      Description
    //  -----   -----  --------------  -------------------------
    //   40     chars                  name - Filename
    //   8      f64    little-endian   mod_time - Seconds since modification
    //   8      u64    little-endian   bytes - File size in bytes
    //
    // ContigRanges Structure (28 bytes)
    //  Bytes   Type   Byte_order      Description
    //  -----   -----  --------------  -------------------------
    //   4      u32    little-endian   bed_idx - Index into BGZ_File
    //   8      u64    little-endian   start_bp - Range start (zero-based)
    //   8      u64    little-endian   end_bp - Range end (zero-based, half-open)
    //   8      u64    little-endian   bgzf_pos - Virtual pos for start of record
    //
    // BGZFiles starts at: 20+(Contigs*44)+(Tiles*4)
    // ContigRanges starts at: 20+(Files*56)+(Contigs*44)+(Tiles*4)
    //
    // Save the ContigIndex to a binary file
    #[allow(dead_code)]
    fn save_index(&self, file_path: &str) -> std::io::Result<()> {
        let fobj = File::create(file_path)?;
        let mut file = io::BufWriter::new(fobj);

        // Write the magic number for this filetype (6-bytes)
        file.write_all(&MAGIC_NUMBER)?;
        // Write the file format version (2-bytes, little-endian)
        file.write_all(&FORMAT_VERSION.to_le_bytes())?;

        // Write tile_size to the file
        file.write_all(&self.tile_size.to_le_bytes())?;

        // Write the number of contigs
        file.write_all(&(self.contigs.len() as u32).to_le_bytes())?;

        // Write the number of bgz files
        file.write_all(&(self.bgz_files.len() as u32).to_le_bytes())?;

        // For each contig write the number of tiles it contains to the file
        for contig in &self.contigs {
            file.write_all(&(contig.contig_tiles.len() as u32).to_le_bytes())?;
        }

        // For each contig/tile write out how many ranges are contained
        for contig in &self.contigs {
            for tile in &contig.contig_tiles {
                file.write_all(&(tile.contig_ranges.len() as u32).to_le_bytes())?;
            }
        }

        for contig in &self.contigs {
            let mut name_vec = vec![0; 40];
            let mut i = 0;
            // TODO: Warn user of name length limitation if name is longer
            //       and make this a program constant tied to the format version
            for byte in contig.name.as_bytes().iter() {
                name_vec[i] = *byte;
                i += 1;
                if i == 40 {
                    break;
                }
            }
            file.write_all(&name_vec)?;
        }

        for bgz_file in &self.bgz_files {
            let mut name_vec = vec![0; 40];
            let mut i = 0;
            // TODO: Warn user of name length limitation if name is longer
            //       and make this a program constant tied to the format version
            // TODO: Consider extending with CRC for file change detection
            for byte in bgz_file.name.as_bytes().iter() {
                name_vec[i] = *byte;
                i += 1;
                if i == 40 {
                    break;
                }
            }
            file.write_all(&name_vec)?;
            file.write_all(&bgz_file.mod_time.to_le_bytes())?;
            file.write_all(&bgz_file.bytes.to_le_bytes())?;
        }

        for contig in &self.contigs {
            for tile in &contig.contig_tiles {
                let mut sorted_ranges = tile.contig_ranges.clone();
                sorted_ranges.sort_by_key(|r| r.start_bp);

                for range in &sorted_ranges {
                    // Write the ContigRange fields
                    file.write_all(&range.bed_idx.to_le_bytes())?;
                    file.write_all(&range.start_bp.to_le_bytes())?;
                    file.write_all(&range.end_bp.to_le_bytes())?;
                    file.write_all(&range.bgzf_pos.to_le_bytes())?;
                }
            }
        }
        Ok(())
    }

    // Save the ContigIndex to a binary file
    #[allow(dead_code)]
    fn save_igd_format(&self, file_path: &str) -> std::io::Result<()> {
        let fobj = File::create(file_path)?;
        let mut file = io::BufWriter::new(fobj);

        // Write tile_size to the file
        file.write_all(&self.tile_size.to_le_bytes())?;

        // Write gType = 1
        let g_type = 1 as i32;
        file.write_all(&g_type.to_le_bytes())?;

        // Write the number of contigs
        file.write_all(&(self.contigs.len() as u32).to_le_bytes())?;

        // For each contig write the number of tiles it contains to the file
        for contig in &self.contigs {
            info!("contig: {} = {}", contig.name, contig.contig_tiles.len());
            file.write_all(&(contig.contig_tiles.len() as u32).to_le_bytes())?;
        }

        // For each contig/tile write out how many ranges are contained
        for contig in &self.contigs {
            for tile in &contig.contig_tiles {
                file.write_all(&(tile.contig_ranges.len() as u32).to_le_bytes())?;
            }
        }

        for contig in &self.contigs {
            let mut name_vec = vec![0; 40];
            let mut i = 0;
            for byte in contig.name.as_bytes().iter() {
                name_vec[i] = *byte;
                i += 1;
                if i == 40 {
                    break;
                }
            }
            file.write_all(&name_vec)?;
        }

        for contig in &self.contigs {
            // Write the number of tiles in the contig
            file.write_all(&(contig.contig_tiles.len() as u32).to_le_bytes())?;

            for tile in &contig.contig_tiles {
                let mut sorted_ranges = tile.contig_ranges.clone();
                //sorted_ranges.sort_by_key(|r| std::cmp::Reverse(r.start_bp));
                sorted_ranges.sort_by_key(|r| r.start_bp);

                let padding = 0 as i32;
                for range in &sorted_ranges {
                    // Write the ContigRange fields
                    file.write_all(&range.bed_idx.to_le_bytes())?;
                    file.write_all(&range.start_bp.to_le_bytes())?;
                    file.write_all(&range.end_bp.to_le_bytes())?;
                    file.write_all(&range.bgzf_pos.to_le_bytes())?;
                    // C-like padding
                    // NOTE: The C compiler will insert padding at the
                    //       end of structures to align the records on
                    //       even word boundaries.  In this case it was
                    //       padding out the 28 byte data to 32 bits.
                    //       There are probably more elegant ways to
                    //       handle this compatibility-mode in Rust.
                    file.write_all(&padding.to_le_bytes())?;
                }
            }
        }
        Ok(())
    }
}

// #[allow(dead_code)]
pub fn prep_idx(
    proj_dir: &String,
    data_type: &String,
) -> Result<(Vec<String>, String, ContigIndex, String), Box<dyn Error>> {
    // Initial instantiation
    let contig_index = ContigIndex {
        tile_size: 16384,
        bgz_files: Vec::new(),
        contig_count: 0,
        contig_lookup: HashMap::new(),
        tile_counts: Vec::new(),
        range_counts: Vec::new(),
        range_data_index: Vec::new(),
        contigs: Vec::new(),
    };

    // TODO: Command line parameter
    // The full directory takes ~4.4 minutes to index
    // The minimal beds take ~57sec to index

    // From the project directory several things can be assumed:
    let index_file = format!("{}/{}_idx.dat", proj_dir, data_type);
    let bgz_dir = format!("{}/{}", proj_dir, data_type);

    let mut filenames: Vec<String> = Vec::new();
    if let Ok(entries) = fs::read_dir(bgz_dir.clone()) {
        for entry in entries {
            if let Ok(entry) = entry {
                if let Some(file_name) = entry.file_name().to_str() {
                    if file_name.ends_with(".bgz") {
                        filenames.push(file_name.to_string());
                    }
                }
            }
        }
    } else {
        error!("Error reading {} directory", bgz_dir);
    }
    Ok((filenames, bgz_dir, contig_index, index_file))
}

#[allow(dead_code)]
pub fn build_idx(
    filenames: &Vec<String>,
    bgz_dir: &String,
    contig_index: &mut ContigIndex,
    index_file: &String,
) {
    let mut fidx = 0;
    for filename in filenames {
        let bgz_file = format!("{}/{}", bgz_dir, filename);
        // Get metadata for the file
        let metadata = fs::metadata(&bgz_file).unwrap();
        // Obtain modification time
        let modification_time = metadata.modified().unwrap();
        // Convert modification time to a more readable format
        let mod_time = modification_time
            .duration_since(SystemTime::UNIX_EPOCH)
            .unwrap();
        // Obtain file size
        let file_size = metadata.len();
        info!(
            "Indexing {} : {} mod_time={:?}, bytes={}",
            fidx, bgz_file, mod_time, file_size
        );
        contig_index.bgz_files.push(BGZFile {
            name: filename.clone(),
            mod_time: mod_time.as_secs_f64(),
            bytes: file_size,
        });

        let mut reader = File::open(bgz_file).map(bgzf::Reader::new).unwrap();
        let mut line = String::new();

        // TODO: validate that the file is in BED format before including it

        let mut virt_pos = u64::from(reader.virtual_position());
        while reader.read_line(&mut line).unwrap() > 0 {
            let fields: Vec<&str> = line.trim_end().split('\t').collect();

            contig_index.add_contig_range(
                fields[0],
                fidx,
                fields[1].parse::<u64>().unwrap(),
                fields[2].parse::<u64>().unwrap(),
                virt_pos,
            );

            // TODO: flush when memory fills.  IGD saves each tile to a file and appends data as it
            // continues.
            // I didn't initially implement this because our use-case doesn't typically challenge the
            // memory of most systems.

            line.clear();
            virt_pos = u64::from(reader.virtual_position());
        }
        fidx += 1;
    }

    let _ = contig_index.save_index(&index_file);
}

#[allow(dead_code)]
pub fn search_idx(
    filenames: &Vec<String>,
    bgz_dir: &String,
    contig_index: &mut ContigIndex,
    index_file: &String,
    q_contig: &String,
    start: u64,
    end: u64,
    family: &Option<String>,
    nrph: bool,
    prod: bool,
) -> Result<Vec<String>, Box<dyn Error>> {
    log::set_logger(&MY_LOGGER).unwrap();
    if prod {
        log::set_max_level(LevelFilter::Warn);
    } else {
        log::set_max_level(LevelFilter::Info);
    }

    debug!("Loading index");
    contig_index.init_search(&index_file);

    // Sanity checking index vs file system
    let mut f_lookup = HashSet::new();
    for filename in filenames {
        f_lookup.insert(filename);
    }
    for ifile in &contig_index.bgz_files {
        if !f_lookup.contains(&ifile.name) {
            warn!("It appears that {} has been deleted from the alignments folder since the index was created!", ifile.name);
        } else {
            f_lookup.remove(&ifile.name);
            let bgz_file = format!("{}/{}", bgz_dir, ifile.name);
            // Get metadata for the file
            let metadata = fs::metadata(bgz_file.clone()).unwrap();
            // Obtain modification time
            let modification_time = metadata.modified().unwrap();
            // Convert modification time to a more readable format
            let mod_time = modification_time
                .duration_since(SystemTime::UNIX_EPOCH)
                .unwrap();
            // Obtain file size
            let file_size = metadata.len();
            if file_size != ifile.bytes {
                warn!("It appears that {} has been modified since the index was created. Byte size difference index={}, file={}", ifile.name, ifile.bytes, file_size);
            } else if mod_time.as_secs_f64() != ifile.mod_time {
                warn!("It appears that {} has been modified since the index was created. Modification time difference index={:?}, file={:?}", ifile.name, ifile.mod_time, mod_time);
            }
        }
    }
    for fsfile in &f_lookup {
        warn!(
            "It appears that {} has been added since the index was created!",
            fsfile
        );
    }

    let mut i_file = File::open(index_file).unwrap();
    debug!("Searching...");
    let results = contig_index.search(&mut i_file, &bgz_dir, &q_contig, start, end, family, nrph);
    return results;
}
