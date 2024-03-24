use crate::derivation::{Literal, Rule};
use anyhow::Result;
use serde::Serialize;
use std::{collections::BTreeMap, fmt::{self, Display, Formatter}, io::BufRead};

#[derive(Debug, Clone, Serialize)]
pub struct Probability {
    lower_bound: f64,
    upper_bound: f64,
}

pub type LiteralProbabilityMap = BTreeMap<Literal, Probability>;
pub type RuleProbabilityMap = BTreeMap<Rule, Probability>;

impl Probability {
    pub const ZERO: Probability = Probability {
        lower_bound: 0.0,
        upper_bound: 0.0,
    };

    pub const ONE: Probability = Probability {
        lower_bound: 1.0,
        upper_bound: 1.0,
    };

    fn new(lower_bound: f64, upper_bound: f64) -> Self {
        Probability {
            lower_bound,
            upper_bound,
        }
    }

    pub fn conjunction(&self, operand: &Probability) -> Probability {
        let lower_bound = 0.0_f64.max(self.lower_bound + operand.lower_bound - 1.0);
        let upper_bound = self.upper_bound.min(operand.upper_bound);
        Probability::new(lower_bound, upper_bound)
    }

    pub fn disjunction(&self, operand: &Probability) -> Probability {
        let lower_bound = self.lower_bound.max(operand.lower_bound);
        let upper_bound = 1.0_f64.min(self.upper_bound + operand.upper_bound);
        Probability::new(lower_bound, upper_bound)
    }

    pub fn multiply(&self, operand: &Probability) -> Probability {
        let lower_bound = self.lower_bound * operand.lower_bound;
        let upper_bound = self.upper_bound * operand.upper_bound;
        Probability::new(lower_bound, upper_bound)
    }

    pub fn try_from_file(path: &str) -> Result<(LiteralProbabilityMap, RuleProbabilityMap)> {
        let file = std::fs::File::open(path)?;
        let reader = std::io::BufReader::new(&file);
        let mut literal_map = LiteralProbabilityMap::new();
        let mut rule_map = RuleProbabilityMap::new();
        reader.lines().for_each(|line| {
            if let Ok(line) = line {
                let parts = line.split_whitespace().collect::<Vec<&str>>();
                if let (Some(&typ), Some(&name), Some(&attributes), Some(&probability)) =
                    (parts.get(0), parts.get(1), parts.get(2), parts.get(3))
                {
                    let attributes_iter = attributes.split(',');
                    let probability = if probability.contains(',') {
                        let probabilities = probability
                            .split(',')
                            .map(|s| s.parse::<f64>().expect("Invalid probability"))
                            .collect::<Vec<f64>>();
                        Probability::new(probabilities[0], probabilities[1])
                    } else {
                        let probability = probability.parse::<f64>().expect("Invalid probability");
                        Probability::new(probability, probability)
                    };

                    match typ {
                        "relation" => {
                            let attributes = attributes_iter
                                .map(|s| s.parse::<u32>().expect("Invalid attribute value"))
                                .collect();
                            let literal = Literal::new(name.to_owned(), attributes);
                            literal_map.insert(literal, probability);
                        }
                        "rule" => {
                            let attributes: Vec<String> =
                                attributes_iter.map(|s| s.to_owned()).collect();
                            let rule = Rule::new(name.to_owned(), attributes);
                            rule_map.insert(rule, probability);
                        }
                        _ => (),
                    }
                }
            }
        });
        Ok((literal_map, rule_map))
    }
}

impl Display for Probability {
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        write!(f, "({}, {})", self.lower_bound, self.upper_bound)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn try_from_file() -> Result<()> {
        let path = "tests/probability.txt";
        let map = Probability::try_from_file(path)?;
        println!("Relations: {:?}\nRules: {:?}", map.0, map.1);
        Ok(())
    }
}
