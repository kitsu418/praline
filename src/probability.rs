use crate::derivation::Literal;
use anyhow::Result;
use csv;
use std::collections::BTreeMap;
use walkdir::WalkDir;

#[derive(Debug)]
struct Probability {
    lower_bound: f32,
    upper_bound: f32,
}

type ProbabilityMap = BTreeMap<Literal, Probability>;

impl Probability {
    fn new(lower_bound: f32, upper_bound: f32) -> Self {
        Probability {
            lower_bound,
            upper_bound,
        }
    }

    fn conjunction(&self, operand: &Probability) -> Probability {
        let lower_bound = 1.0
            - operand.upper_bound.min(1.0 - self.lower_bound)
            - self.upper_bound.min(1.0 - operand.lower_bound)
            - (1.0 - self.lower_bound).min(1.0 - operand.lower_bound);
        let upper_bound = self.upper_bound.min(operand.upper_bound);
        Probability::new(lower_bound, upper_bound)
    }

    fn disjunction(&self, operand: &Probability) -> Probability {
        let lower_bound = self.lower_bound.max(operand.lower_bound);
        let upper_bound = 1.0_f32.min(self.upper_bound + operand.upper_bound);
        Probability::new(lower_bound, upper_bound)
    }

    fn try_from_csv(path: &str) -> Result<ProbabilityMap> {
        let mut reader = csv::Reader::from_path(path)?;
        let map: ProbabilityMap = reader
            .records()
            .filter_map(|result| {
                result.ok().map(|record| {
                    // println!("{:?}", record);
                    let relation_name = record.get(0).expect("Missing relation name").to_owned();
                    let attributes = record
                        .iter()
                        .skip(1)
                        .take(record.len() - 3)
                        .map(|s| s.parse::<u32>().expect("Invalid attribute value"))
                        .collect();
                    let literal = Literal::new(relation_name, attributes);

                    let lower_bound = record
                        .get(record.len() - 2)
                        .expect("Missing lower bound")
                        .parse::<f32>()
                        .expect("Invalid lower bound");
                    let upper_bound = record
                        .get(record.len() - 1)
                        .expect("Missing upper bound")
                        .parse::<f32>()
                        .expect("Invalid upper bound");
                    let probability = Probability::new(lower_bound, upper_bound);

                    (literal, probability)
                })
            })
            .collect();
        Ok(map)
    }

    fn try_from_dir(path: &str) -> Result<ProbabilityMap> {
        let map = WalkDir::new(path)
            .into_iter()
            .fold(ProbabilityMap::new(), |mut acc, entry| {
                if let Ok(entry) = entry {
                    if entry.file_type().is_file() {
                        if entry.path().extension().unwrap() == "csv" {
                            let path = entry.path().to_str().expect("Invalid path");
                            let map = Probability::try_from_csv(path).expect("Invalid csv");
                            acc.extend(map);
                        }
                    }
                }
                acc
            });
        Ok(map)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_try_from_csv_path() -> Result<()> {
        let path = "tests/probability/edge.csv";
        let map = Probability::try_from_csv(path)?;
        println!("{:?}", map);
        Ok(())
    }

    #[test]
    fn test_try_from_dir_path() -> Result<()> {
        let path = "tests/probability";
        let map = Probability::try_from_dir(path)?;
        println!("{:?}", map);
        Ok(())
    }
}
