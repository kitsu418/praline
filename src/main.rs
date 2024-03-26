use anyhow::Result;
use clap::Parser;
use problog::analysis;
use problog::derivation::Derivation;
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
    let probability = Probability::load(&args.probability)?;
    let mut analysis = analysis::Analysis::new(derivation, probability.0, probability.1);
    analysis.calculate_probability();
    analysis.dump(args.output)?;
    Ok(())
}
