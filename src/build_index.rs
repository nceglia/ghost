//! BuildIndex stage code

use serde::{Serialize, Deserialize};
use martian::prelude::*;
use martian_derive::*;
use bio::io::fasta;
use debruijn_mapping::{
    build_index::build_index,
    mappability::analyze_graph,
};
use debruijn_mapping::{config, utils};
use std::path::PathBuf;

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct BuildIndexStageInputs {
    ref_fasta: PathBuf,
}

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct BuildIndexStageOutputs {
    index: PathBuf,
}

pub struct BuildIndex;

#[make_mro(mem_gb = 8, threads = 20)]
impl MartianMain for BuildIndex {
    type StageInputs = BuildIndexStageInputs;
    type StageOutputs = BuildIndexStageOutputs;
    fn main(
        &self,
        args: Self::StageInputs,
        rover: MartianRover,
    ) -> Result<Self::StageOutputs, Error> {
        let index_file: PathBuf = rover.make_path("transcriptome.idx");
        println!("Reading Reference.");
        let reference = fasta::Reader::from_file(args.ref_fasta)?;
        let (seqs, tx_names, tx_gene_map) = utils::read_transcripts(reference)?;
        let index = build_index::<config::KmerType>(&seqs, &tx_names, &tx_gene_map, 20)?;
        println!("Writing index to disk");
        utils::write_obj(&index, &index_file)?;
        println!("Finished writing index!");
        let index = debruijn_mapping::utils::read_obj(String::from("transcriptome.idx"))?;
        let records = analyze_graph::<config::KmerType>(&index)?;
        println!("{} transcripts total", records.len());
        Ok(BuildIndexStageOutputs {
            index: index_file,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn run_stage() {
        let args = BuildIndexStageInputs {
            ref_fasta: PathBuf::from("/work/shah/reference/transcriptomes/fastas/gencode_sm.fa"),
        };
        let stage = BuildIndex;
        let res = stage.test_run_tmpdir(args).unwrap();
        assert_eq!(res.index.file_name().is_some(), PathBuf::from("transcriptome.idx").file_name().is_some());
    }
}
