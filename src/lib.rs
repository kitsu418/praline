pub mod derivation;
pub mod probability;
pub mod analysis;

use anyhow::Result;

trait Load {
    fn load<F>(path: &str) -> Result<Self> where Self: Sized;
}