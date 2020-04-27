//! CorrectBarcodes stage code

use serde::{Serialize, Deserialize};
use martian::prelude::*;
use martian_derive::*;
use debruijn::dna_string::DnaString;
use failure::Error;
use bio::io::fastq;
use std::{
   fs::File,
   io::{prelude::*, BufReader},
   path::Path,
};
use std::path::PathBuf;
use rayon::prelude::*;
use std::collections::HashSet;

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct CorrectBarcodesStageInputs {
    fastq_r1: PathBuf,
    whitelist: PathBuf,
}

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct CorrectBarcodesStageOutputs {
    valid_barcodes: Vec<String>,
}

pub struct CorrectBarcodes;

fn read_barcodes(filename: impl AsRef<Path>) -> HashSet<String> {
    let file = File::open(filename).expect("no such file");
    let buffer = BufReader::new(file);
    let barcodes: Vec<String> = buffer.lines().map(|l| l.expect("Could not parse line")).collect();
    let mut barcode_set = HashSet::new();
    for barcode in barcodes {
        barcode_set.insert(barcode);
    }
    barcode_set
 }
 
 fn mismatches(sequence : &String, reference : &String) -> u32 {
    let mut mismatch : u32 = 0;
    for (known_nt, nt) in reference.chars().zip(sequence.chars()) {
        if known_nt != nt {
            mismatch = mismatch + 1;
            if mismatch > 1 {
                break;
            }
        }
    }
    mismatch
 }
 
 fn mismatches_par(sequence: &String, references : HashSet<String>) -> u32 {
    references.par_iter().map(|x| mismatches(sequence, &x)).min().unwrap()
 }
 
 fn parse_barcode(record: &bio::io::fastq::Record) -> String {
    let _seqid = record.id().to_owned();
    let sequence = DnaString::from_acgt_bytes_hashn(record.seq(), record.id().as_bytes());
    let barcode = sequence.slice(0, 16);
    barcode.to_string()
 }

 #[make_mro(mem_gb = 4, threads = 12)]
impl MartianMain for CorrectBarcodes {
    type StageInputs = CorrectBarcodesStageInputs;
    type StageOutputs = CorrectBarcodesStageOutputs;
    fn main(
        &self,
        args: Self::StageInputs,
        _rover: MartianRover,
    ) -> Result<Self::StageOutputs, Error> {
        println!("Reading R1 Fastq...");
        let records = fastq::Reader::from_file(&args.fastq_r1)?;
        println!("Finished Reading R1.\n");

        let barcode_set = read_barcodes(args.whitelist);
        println!("Finished Hashing Barcode Set.\n");

        let mut read_records = Vec::new();

        for record_ref in records.records() {
            let record = record_ref?;
            read_records.push(record.to_owned());
        }

        println!("\nParsing Sequence Barcodes..");
        let barcodes = read_records.par_iter().map(|x| parse_barcode(x));
        println!("Finished.");

        println!("\nFinding Valid Barcode Set...");
        let mut valid_barcodes : Vec<String> = barcodes.clone().filter(|x| barcode_set.contains(x)).collect();
        valid_barcodes.par_sort_unstable();
        valid_barcodes.dedup();
        println!("Perfectly Matched Barcodes: {:?}", valid_barcodes.len());

        println!("\nFinding Invalid Barcode Set...");
        let mut invalid_barcodes : Vec<String> = barcodes.clone().filter(|x| !barcode_set.contains(x)).collect();
        invalid_barcodes.par_sort_unstable();
        invalid_barcodes.dedup();
        println!("Invalid Barcodes: {:?}", invalid_barcodes.len());

        let corrected_barcodes : Vec<String> = barcodes.clone().filter(|x| mismatches_par(x, barcode_set.to_owned()) <= 1u32).collect();

        for barcode in corrected_barcodes {
            if !valid_barcodes.contains(&barcode) {
                valid_barcodes.push(barcode)
            }
        }
        
        println!("\n{:?} - Valid Barcodes", valid_barcodes.len());
        Ok(CorrectBarcodesStageOutputs {
            valid_barcodes: valid_barcodes,
        })
    }
}

#[cfg(test)]
mod tests {     
    use super::*;
    #[test]
    fn run_stage() {
        let args = CorrectBarcodesStageInputs {
            whitelist: PathBuf::from("/work/shah/ceglian/reaper/scdiff/test/3M-february-2018.txt"),
            fastq_r1: PathBuf::from("/work/shah/ceglian/reaper/ghost/data/small_r1.fastq"),
        };
        let stage = CorrectBarcodes;
        let res = stage.test_run_tmpdir(args).unwrap();
        assert_eq!(res.valid_barcodes.len(), 344);
    }
}

