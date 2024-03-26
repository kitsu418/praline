use std::fmt::Display;

use anyhow::Result;
use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize, Debug, Ord, PartialEq, PartialOrd, Eq, Clone)]
pub struct Relation {
    pub name: String,
    attributes: Vec<u32>,
}

type Body = Vec<Relation>;

impl Relation {
    pub fn new(name: String, attributes: Vec<u32>) -> Self {
        Relation {
            name,
            attributes,
        }
    }
}

#[derive(Debug, PartialEq, PartialOrd, Eq, Ord, Clone)]
pub struct Rule {
    pub head: String,
    pub body: Vec<String>,
}

impl Rule {
    pub fn new(head: String, body: Vec<String>) -> Self {
        Rule { head, body }
    }
}

#[derive(Serialize, Deserialize, Debug, Clone, PartialEq, PartialOrd, Eq, Ord)]
pub struct Derivation {
    pub head: Relation,
    pub bodies: Vec<Body>,
}

impl Derivation {
    pub fn try_from_json(path: &str) -> Result<Vec<Self>> {
        let json_file = std::fs::File::open(path)?;
        let json_reader = std::io::BufReader::new(json_file);
        let derivations: Vec<Derivation> = serde_json::from_reader(json_reader)?;
        Ok(derivations)
    }
}

impl Display for Relation {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}(", self.name)?;
        for (i, attribute) in self.attributes.iter().enumerate() {
            write!(f, "{}", attribute)?;
            if i < self.attributes.len() - 1 {
                write!(f, ", ")?;
            }
        }
        write!(f, ")")?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_dump() -> Result<()> {
        let path = "tests/derivation.json";
        let derivations = Derivation::try_from_json(path)?;
        println!("{:?}", derivations);
        Ok(())
    }
}
