
//! Martian-rust adapter ghost

use serde::Deserialize;
use martian::prelude::*;
use docopt::Docopt;

//Stages
mod build_index;
mod map_reads;
mod correct_barcodes;
mod count_umi;

const USAGE: &'static str = "
Martian adapter for ghost executable

Usage:
  ghost martian <adapter>...
  ghost mro [--file=<filename>] [--rewrite]
  ghost --help

Options:
  --help              Show this screen.
  --file=<filename>   Output filename for the mro.
  --rewrite           Whether to rewrite the file if it exists.
";

#[derive(Debug, Deserialize)]
struct Args {
    // Martian interface
    cmd_martian: bool,
    arg_adapter: Vec<String>,
    // Mro generation
    cmd_mro: bool,
    flag_file: Option<String>,
    flag_rewrite: bool,
}

fn main() -> Result<(), Error> {

    let args: Args = Docopt::new(USAGE)
        .and_then(|d| d.deserialize())
        .unwrap_or_else(|e| e.exit());

    let (stage_registry, mro_registry) = martian_stages![
        build_index::BuildIndex,
        map_reads::MapReads,
        correct_barcodes::CorrectBarcodes,
        count_umi::CountUmi,
    ];

    if args.cmd_martian {
        // Setup the martian adapter
        let runner = MartianAdapter::new(stage_registry);
        // If you want explicit control over the log level, use:
        // let runner = runner.log_level();
        // run the stage
        let retcode = runner.run(args.arg_adapter);
        // return from the process
        std::process::exit(retcode);
    } else if args.cmd_mro {
        // Create the mro for all the stages in this adapter
        martian_make_mro(args.flag_file, args.flag_rewrite, mro_registry)?;
    } else {
        // If you need custom commands, implement them here
        unimplemented!()
    }
    
    Ok(())
}
