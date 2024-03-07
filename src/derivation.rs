use anyhow::Result;
use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize, Debug, Ord, PartialEq, PartialOrd, Eq, Clone)]
pub struct Literal {
    relation_name: String,
    attributes: Vec<u32>,
}

type Conjunction = Vec<Literal>;
type Disjunction = Vec<Conjunction>;

impl Literal {
    pub fn new(relation_name: String, attributes: Vec<u32>) -> Self {
        Literal {
            relation_name,
            attributes,
        }
    }
}

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct Derivation {
    pub parent: Literal,
    pub children: Disjunction,
}

impl Derivation {
    pub fn try_from_json(path: &str) -> Result<Vec<Self>> {
        let json_file = std::fs::File::open(path)?;
        let json_reader = std::io::BufReader::new(json_file);
        let derivations: Vec<Derivation> = serde_json::from_reader(json_reader)?;
        Ok(derivations)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_try_from_json_path() -> Result<()> {
        let path = "tests/derivation.json";
        let derivations = Derivation::try_from_json(path)?;
        println!("{:?}", derivations);
        Ok(())
    }
}
