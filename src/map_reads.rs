//! MapReads stage code

use serde::{Serialize, Deserialize};

// The prelude brings the following items in scope:
// - Traits: MartianMain, MartianStage, RawMartianStage, MartianFileType, MartianMakePath
// - Struct/Enum: MartianAdapter, MartianRover, Resource, StageDef, MartianVoid,
//                Error (from failure crate), LevelFilter (from log crate)
// - Macros: martian_stages!
// - Functions: martian_make_mro
use martian::prelude::*;
use std::collections::HashMap;
use martian_derive::*;
use bio::io::fastq;
use std::path::PathBuf;

use debruijn_mapping::{config, utils};

use std::fmt::Debug;
use std::io::{self, Write};
use std::sync::{mpsc, Arc, Mutex};

use crossbeam_utils::thread::scope;
use debruijn::dna_string::DnaString;
use failure::Error;

use config::READ_COVERAGE_THRESHOLD;
use debruijn_mapping::pseudoaligner::Pseudoaligner;

// NOTE: The following two structs will serve as the associated type for the
// trait. The struct fields need to be owned and are limited to
// - Basic int/float/bool/String types, PathBuf, Vec, Option, HashMap, HashSet
// - Structs/Enums implementing "AsMartianPrimaryType" (You can use #[derive(MartianType)])
// - Filetype (see the note below, representing as a filetype in mro)

// If you want to declare a new filetype use the `martian_filetype!` macro:
// martian_filetype!(Lz4File, "lz4");

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct MapReadsStageInputs {
    index: PathBuf,
    fastq_r2: PathBuf,
}

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct MapReadsStageOutputs {
    eq_classes: HashMap::<std::string::String,Vec<u32>>,
    coverage: HashMap::<std::string::String,usize>,
}

pub struct MapReads;

#[make_mro(mem_gb = 4, threads = 12)]
impl MartianMain for MapReads {
    type StageInputs = MapReadsStageInputs;
    type StageOutputs = MapReadsStageOutputs;
    fn main(
        &self,
        args: Self::StageInputs,
        _rover: MartianRover,
    ) -> Result<Self::StageOutputs, Error> {
        let num_threads = 12;
        println!("Reading index from disk");
        let index : Pseudoaligner<debruijn::kmer::Kmer20> = utils::read_obj(args.index)?;
        println!("Mapping reads from fastq");
        let reads = fastq::Reader::from_file(args.fastq_r2)?;
        let mut mapped = Vec::<(bool, std::string::String, std::vec::Vec<u32>, usize)>::new();
        println!("Starting Multi-threaded Mapping");
        let (tx, rx) = mpsc::sync_channel(num_threads);
        let atomic_reader = Arc::new(Mutex::new(reads.records()));
        println!("Spawning {} threads for Mapping.\n", num_threads);
        scope(|scope| {
            for _ in 0..num_threads {
                let tx = tx.clone();
                let reader = Arc::clone(&atomic_reader);
                let index = &index;
                scope.spawn(move |_| {
                    loop {
                        match utils::get_next_record(&reader) {
                            Some(result_record) => {
                                let record = match result_record {
                                    Ok(record) => record,   
                                    Err(err) => panic!("Error {:?} in reading fastq", err),
                                };
                                let dna_string = std::str::from_utf8(record.seq()).unwrap();
                                let seq = DnaString::from_dna_string(dna_string);
                                let read_data = index.map_read(&seq);
                                let wrapped_read_data = match read_data {
                                    Some((eq_class, coverage)) => {
                                        if coverage >= READ_COVERAGE_THRESHOLD && eq_class.is_empty() {
                                            Some((true, record.id().to_owned(), eq_class, coverage))
                                        } else {
                                            Some((false, record.id().to_owned(), eq_class, coverage))
                                        }
                                    }
                                    None => Some((false, record.id().to_owned(), Vec::new(), 0usize)),
                                };
                                tx.send(wrapped_read_data).expect("Could not send data!");
                            }
                            None => {
                                tx.send(None).expect("Could not send data!");
                                break;
                            }
                        };
                    }
                });
            }
            let mut read_counter: usize = 0;
            let mut mapped_read_counter: usize = 0;
            let mut dead_thread_count = 0;
            for eq_class in rx.iter() {
                match eq_class {
                    None => {
                        dead_thread_count += 1;
                        if dead_thread_count == num_threads {
                            drop(tx);
                            for eq_class in rx.iter() {
                                eq_class.map_or((), |eq_class| eprintln!("{:?}", eq_class));
                            }
                            break;
                        }
                    }
                    Some(read_data) => {
                        if read_data.0 {
                            mapped_read_counter += 1;
                        }
                        read_counter += 1;
                        if read_counter % 1000 == 0 {
                            let frac_mapped = mapped_read_counter as f32 * 100.0 / read_counter as f32;
                            eprint!(
                                "\rDone Mapping {} reads w/ Rate: {}",
                                read_counter, frac_mapped
                            );
                            io::stderr().flush().expect("Could not flush stdout");
                        }
                        mapped.push(read_data);
                    } 
                }
            } 
        })
        .unwrap();
        let mut eq_classes = HashMap::new();
        let mut coverage = HashMap::new();
        for read_data in mapped.iter() {
            if read_data.0 {
                println!("{}",read_data.1);
                eq_classes.insert(read_data.1.clone(),read_data.2.clone());
                coverage.insert(read_data.1.clone(),read_data.3.clone());
            }
        }
        Ok(MapReadsStageOutputs {
            eq_classes: eq_classes,
            coverage: coverage,
        })
    }
}

#[cfg(test)]
mod tests {     
    use super::*;
    #[test]
    fn run_stage() {
        let args = MapReadsStageInputs {
            index: PathBuf::from("data/transcriptome.idx"),
            fastq: PathBuf::from("/work/shah/ceglian/reaper/ghost/data/small_r2.fastq"),
        };
        let stage = MapReads;
        let res = stage.test_run_tmpdir(args).unwrap();
        assert_eq!(res.eq_classes.len(), 131);
        assert_eq!(res.coverage.len(), 131);
    }
}

