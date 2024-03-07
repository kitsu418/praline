use std::fs::File;
use std::io::{self, BufWriter, Write};

use anyhow::Result;
use clap::Parser;
use problog::analysis;
use problog::derivation::{Derivation, Literal};
use problog::probability::Probability;

#[derive(Parser)]
#[command(version, about, long_about)]
struct Args {
    /// path to the dumped derivation tree json file
    #[arg(short, long)]
    derivation: String,
    /// path to the probability csv file or the directory containing the probability csv files
    #[arg(short, long)]
    probability: String,
    #[arg(short, long)]
    output: Option<String>,
}

fn main() -> Result<()> {
    let args = Args::parse();
    let derivation = Derivation::try_from_json(&args.derivation)?;
    let probability = Probability::try_from_dir(&args.probability)?;
    let mut analysis = analysis::Analysis::new(derivation, probability);
    analysis.calculate_probability();

    let writer: Box<dyn Write> = match args.output {
        Some(output) => {
            let file = File::create(output).expect("Failed to create output file");
            Box::new(BufWriter::new(file))
        }
        None => Box::new(BufWriter::new(io::stdout().lock())),
    };

    let probability_pairs: Vec<(Literal, Probability)> =
        analysis.probability_map.clone().into_iter().collect();
    serde_json::to_writer(writer, &probability_pairs)
        .expect("Failed to dump result to the output file.");
    Ok(())
}
