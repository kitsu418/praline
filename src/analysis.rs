use std::collections::VecDeque;

use crate::derivation::Derivation;
use crate::probability::{Probability, ProbabilityMap};

pub struct Analysis {
    derivations: Vec<Derivation>,
    pub probability_map: ProbabilityMap,
}

impl Analysis {
    pub fn new(derivations: Vec<Derivation>, probability_map: ProbabilityMap) -> Self {
        Analysis {
            derivations,
            probability_map,
        }
    }

    pub fn calculate_probability(&mut self) {
        let mut worklist: VecDeque<&Derivation> = self
            .derivations
            .iter()
            .filter(|&derivation| !self.probability_map.contains_key(&derivation.parent))
            .collect();
        while let Some(derivation) = worklist.pop_front() {
            if let Ok(probability) =
                derivation
                    .children
                    .iter()
                    .try_fold(Probability::ZERO, |acc, conjunction| {
                        conjunction
                            .iter()
                            .try_fold(Probability::ONE, |acc, literal| {
                                self.probability_map
                                    .get(literal)
                                    .map_or(Err(()), |p| Ok(acc.conjunction(p)))
                            })
                            .map_or(Err(()), |p| Ok(acc.disjunction(&p)))
                    })
            {
                self.probability_map
                    .insert(derivation.parent.clone(), probability);
            } else {
                worklist.push_back(derivation);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use anyhow::Result;

    #[test]
    fn test_calculate_probability() -> Result<()> {
        let derivations = Derivation::try_from_json("tests/derivation.json").unwrap();
        let probability_map = Probability::try_from_csv("tests/probability/edge.csv").unwrap();
        let mut analysis = Analysis::new(derivations, probability_map);
        analysis.calculate_probability();
        println!("{:?}", analysis.probability_map);
        Ok(())
    }
}
