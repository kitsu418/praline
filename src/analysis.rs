use crate::derivation::{Derivation, Relation, Rule};
use crate::probability::{Probability, RelationProbabilityMap, RuleProbabilityMap};
use anyhow::Result;
use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufWriter, Write};

pub struct Analysis {
    derivations: Vec<Derivation>,
    pub relation_probability_map: RelationProbabilityMap,
    rule_probability_map: RuleProbabilityMap,
}

#[derive(Debug)]
enum State {
    Known(Probability),
    Unknown(Vec<usize>),
}

impl Analysis {
    pub fn new(
        derivations: Vec<Derivation>,
        relation_probability_map: RelationProbabilityMap,
        rule_probability_map: RuleProbabilityMap,
    ) -> Self {
        Analysis {
            derivations,
            relation_probability_map,
            rule_probability_map,
        }
    }

    fn get_rule_probability(&self, head: &Relation, body: &[Relation]) -> Option<&Probability> {
        let rule = Rule::new(
            head.name.clone(),
            body.iter().map(|r| r.name.clone()).collect(),
        );
        self.rule_probability_map.get(&rule)
    }

    pub fn calculate_probability(&mut self) {
        let mut state_map: BTreeMap<(&Relation, &Vec<Relation>), State> = BTreeMap::new();
        let mut fixed = false;
        let worklist: Vec<&Derivation> = self
            .derivations
            .iter()
            .filter_map(|d| {
                if self.relation_probability_map.contains_key(&d.head) {
                    None
                } else {
                    Some(d)
                }
            })
            .collect();
        while !fixed {
            fixed = true;
            for derivation in &worklist {
                if self.relation_probability_map.contains_key(&derivation.head) {
                    continue;
                }
                let mut probability = Some(Probability::ZERO);
                for body in &derivation.bodies {
                    if let Some(State::Known(child_probability)) =
                        state_map.get(&(&derivation.head, body))
                    {
                        probability = probability.map(|p| p.disjunction(child_probability));
                        continue;
                    }

                    let mut body_probability = Probability::ONE;
                    let mut unknown = Vec::<usize>::new();

                    for (id, relation) in body.iter().enumerate() {
                        if let Some(p) = self.relation_probability_map.get(relation) {
                            body_probability = body_probability.conjunction(p);
                        } else {
                            unknown.push(id);
                        }
                    }

                    if unknown.is_empty() {
                        if let Some(rule_probability) =
                            self.get_rule_probability(&derivation.head, body)
                        {
                            body_probability = body_probability.multiply(rule_probability)
                        }
                        if probability.is_some() {
                            probability = probability.map(|p| p.disjunction(&body_probability));
                        }
                        state_map.insert((&derivation.head, body), State::Known(body_probability));
                        fixed = false;
                    } else {
                        probability = None;
                        if let Some(State::Unknown(previous)) =
                            state_map.get(&(&derivation.head, body))
                        {
                            if previous.eq(&unknown) {
                                continue;
                            }
                        }
                        state_map.insert((&derivation.head, body), State::Unknown(unknown));
                        fixed = false;
                    }
                }
                if let Some(probability) = probability {
                    self.relation_probability_map
                        .insert(derivation.head.clone(), probability);
                }
            }
        }
    }

    pub fn dump(&self, output: Option<String>) -> Result<()> {
        let mut writer: Box<dyn Write> = match output {
            Some(output) => {
                let file = File::create(output).expect("Failed to create output file");
                Box::new(BufWriter::new(file))
            }
            None => Box::new(BufWriter::new(std::io::stdout().lock())),
        };
        for (relation, probability) in &self.relation_probability_map {
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
        let (relation_probability_map, rule_probability_map) =
            Probability::load("tests/probability.txt").unwrap();
        let mut analysis =
            Analysis::new(derivations, relation_probability_map, rule_probability_map);
        analysis.calculate_probability();
        analysis.dump(None)?;
        Ok(())
    }
}
