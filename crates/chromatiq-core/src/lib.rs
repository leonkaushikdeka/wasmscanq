//! Chromatiq Core Library
//!
//! A high-performance genomic data processing library that provides
//! binary parsing, pileup generation, and coverage calculation algorithms.
//!
//! # Example
//!
//! ```ignore
//! use chromatiq_core::{BamReader, Coverage};
//!
//! let reader = BamReader::new("sample.bam")?;
//! let coverage = reader.calculate_coverage(0, 1000000)?;
//! ```

#![cfg_attr(not(feature = "std"), no_std)]
#![warn(missing_docs)]

extern crate alloc;

mod bam;
mod coverage;
mod pileup;
mod reference;
mod errors;

pub use bam::{BamReader, BamRecord, ReadFlags};
pub use coverage::Coverage;
pub use pileup::{Pileup, PileupColumn};
pub use reference::{ReferenceSequence, GenomeAssembly};
pub use errors::{CoreError, CoreResult};
