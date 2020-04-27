
//! CountUmi stage code

use serde::{Serialize, Deserialize};
use martian::prelude::*;
use martian_derive::*;
use debruijn::dna_string::DnaString;
use failure::Error;
use bio::io::fastq;
use std::path::PathBuf;
use rayon::prelude::*;

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct CountUmiStageInputs {
    fastq_r1: PathBuf,
}

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct CountUmiStageOutputs {
    umis: Vec<String>,
}

fn parse_umi(record: &bio::io::fastq::Record) -> String {
    let _seqid = record.id().to_owned();
    let sequence = DnaString::from_acgt_bytes_hashn(record.seq(), record.id().as_bytes());
    let umi = sequence.slice(16, 26);
    umi.to_string()
 }

// This is our stage struct
pub struct CountUmi;

#[make_mro(mem_gb = 4, threads = 4)]
impl MartianMain for CountUmi {
    type StageInputs = CountUmiStageInputs;
    type StageOutputs = CountUmiStageOutputs;
    fn main(
        &self,
        args: Self::StageInputs,
        _rover: MartianRover,
    ) -> Result<Self::StageOutputs, Error> {
        println!("Reading R1.\n");
        let records = fastq::Reader::from_file(&args.fastq_r1)?;
        println!("Finished Reading R1.\n");

        let mut read_records = Vec::new();

        for record_ref in records.records() {
            let record = record_ref?;
            read_records.push(record.to_owned());
        }

        println!("\nParsing Sequence UMIs...");
        let umis = read_records.par_iter().map(|x| parse_umi(x)).collect();
        println!("Finished.");
        Ok(CountUmiStageOutputs {
            umis: umis,
        })
    }
}

#[cfg(test)]
mod tests {     
    use super::*;
    #[test]
    fn run_stage() {
        let args = CountUmiStageInputs {
            fastq_r1: PathBuf::from("/work/shah/ceglian/reaper/ghost/data/small_r1.fastq"),
        };
        let stage = CountUmi;
        let res = stage.test_run_tmpdir(args).unwrap();
        assert_eq!(res.umis.len(), 750);
    }
}
