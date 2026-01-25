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
mod errors;
mod indexed_bam;
mod pileup;
mod reference;
mod vcf;

pub use bam::{BamReader, BamRecord, CigarOp, CigarOperation, ReadFlags, ReferenceInfo};
pub use coverage::{
    calculate_coverage_downsampled, calculate_coverage_simple, Coverage, CoverageCalculator,
};
pub use errors::{CoreError, CoreResult};
pub use indexed_bam::{BaiIndex, IndexedBamReader, RegionIterator, RegionStats};
pub use pileup::{
    find_variants, generate_pileup, Pileup, PileupColumn, PileupConfig, PileupGenerator,
    VariantCandidate,
};
pub use reference::{GenomeAssembly, GenomicRegion, ReferenceCollection, ReferenceSequence};
pub use vcf::{VariantType, VcfHeader, VcfReader, VcfRecord, VcfStats};
