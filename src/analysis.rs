use crate::derivation::{Derivation, Literal, Rule};
use crate::probability::{LiteralProbabilityMap, Probability, RuleProbabilityMap};
use anyhow::Result;
use std::collections::{BTreeMap, VecDeque};
use std::fs::File;
use std::io::{BufWriter, Write};

pub struct Analysis {
    derivations: Vec<Derivation>,
    pub literal_probability_map: LiteralProbabilityMap,
    rule_probability_map: RuleProbabilityMap,
}

enum State {
    Known(Probability),
    Unknown(Vec<usize>),
}

impl Analysis {
    pub fn new(
        derivations: Vec<Derivation>,
        literal_probability_map: LiteralProbabilityMap,
        rule_probability_map: RuleProbabilityMap,
    ) -> Self {
        Analysis {
            derivations,
            literal_probability_map,
            rule_probability_map,
        }
    }

    fn get_rule_probability(&self, head: &Literal, body: &[Literal]) -> Option<&Probability> {
        let rule = Rule::new(
            head.relation_name.clone(),
            body.iter().map(|l| l.relation_name.clone()).collect(),
        );
        self.rule_probability_map.get(&rule)
    }

    pub fn calculate_probability(&mut self) {
        let mut worklist: VecDeque<&Derivation> = self
            .derivations
            .iter()
            .filter(|&derivation| {
                !self
                    .literal_probability_map
                    .contains_key(&derivation.parent)
            })
            .collect();
        let mut state_map: BTreeMap<(&Literal, &Vec<Literal>), State> = BTreeMap::new();
        while let Some(derivation) = worklist.pop_front() {
            let mut probability = Some(Probability::ZERO);
            let mut fixed = true;
            for conjunction in &derivation.children {
                if let Some(State::Known(child_probability)) =
                    state_map.get(&(&derivation.parent, conjunction))
                {
                    probability = probability.map(|p| p.disjunction(child_probability));
                }

                let mut child_probability = Probability::ONE;
                let mut unknown = Vec::<usize>::new();

                for (id, literal) in conjunction.iter().enumerate() {
                    if let Some(p) = self.literal_probability_map.get(literal) {
                        child_probability = child_probability.conjunction(p);
                    } else {
                        unknown.push(id);
                    }
                }
                if unknown.is_empty() {
                    child_probability = if let Some(rule_probability) =
                        self.get_rule_probability(&derivation.parent, conjunction)
                    {
                        child_probability.multiply(rule_probability)
                    } else {
                        child_probability
                    };
                    if probability.is_some() {
                        probability = probability.map(|p| p.disjunction(&child_probability));
                    }
                    state_map.insert(
                        (&derivation.parent, conjunction),
                        State::Known(child_probability),
                    );
                } else {
                    probability = None;
                    if let Some(State::Unknown(previous)) =
                        state_map.get(&(&derivation.parent, conjunction))
                    {
                        if !previous.eq(&unknown) {
                            fixed = false;
                        }
                    } else {
                        fixed = false;
                    }
                    state_map.insert((&derivation.parent, conjunction), State::Unknown(unknown));
                }
            }
            if let Some(probability) = probability {
                self.literal_probability_map
                    .insert(derivation.parent.clone(), probability);
            } else if !fixed {
                worklist.push_back(derivation);
            }
        }
        for ((head, body), state) in state_map {
            if let State::Unknown(_) = state {
                println!("{:?} {:?}", head, body);
            }
        }
    }

    pub fn try_dump_probability_map(&self, output: Option<String>) -> Result<()> {
        let mut writer: Box<dyn Write> = match output {
            Some(output) => {
                let file = File::create(output).expect("Failed to create output file");
                Box::new(BufWriter::new(file))
            }
            None => Box::new(BufWriter::new(std::io::stdout().lock())),
        };
        for (relation, probability) in &self.literal_probability_map {
            writer
                .write_all(format!("{}: {}\n", relation, probability).as_bytes())
                .expect("Failed to write to the output file.");
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use anyhow::Result;

    #[test]
    fn test_calculate_probability() -> Result<()> {
        let derivations = Derivation::try_from_json("tests/derivation.json").unwrap();
        let (literal_probability_map, rule_probability_map) =
            Probability::try_from_file("tests/probability.txt").unwrap();
        let mut analysis =
            Analysis::new(derivations, literal_probability_map, rule_probability_map);
        analysis.calculate_probability();
        println!("{:?}", analysis.try_dump_probability_map(None)?);
        Ok(())
    }
}
