//! Reference sequence handling
//!
//! This module provides types and utilities for managing reference genome
//! sequences and their metadata.

use alloc::string::String;
use alloc::vec::Vec;

/// Supported genome assemblies
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GenomeAssembly {
    /// GRCh38 (Human Reference Genome 38)
    GRCh38,
    /// GRCh37 (Human Reference Genome 37)
    GRCh37,
    /// Mouse reference genome
    GRCm38,
    /// Custom assembly
    Custom,
}

impl GenomeAssembly {
    /// Get the display name for the assembly
    pub fn display_name(&self) -> &'static str {
        match self {
            GenomeAssembly::GRCh38 => "GRCh38",
            GenomeAssembly::GRCh37 => "GRCh37",
            GenomeAssembly::GRCm38 => "GRCm38",
            GenomeAssembly::Custom => "Custom",
        }
    }
}

/// Reference sequence information
#[derive(Debug, Clone)]
pub struct ReferenceSequence {
    /// Accession number (e.g., "NC_000001.11")
    pub accession: String,
    /// Display name (e.g., "chr1")
    pub name: String,
    /// Alternative names
    pub aliases: Vec<String>,
    /// Sequence length in bases
    pub length: u32,
    /// MD5 checksum of the sequence
    pub md5: Option<String>,
    /// Genome assembly
    pub assembly: GenomeAssembly,
    /// Species (e.g., "Homo sapiens")
    pub species: String,
    /// Is this a chromosome (vs contig/scaffold)
    pub is_chromosome: bool,
    /// Chromosome name (for mitochondria, etc.)
    pub chromosome: Option<String>,
}

impl ReferenceSequence {
    /// Create a new reference sequence
    pub fn new(accession: String, name: String, length: u32, species: String) -> Self {
        let is_chromosome = name.starts_with("chr") && name[3..].parse::<u32>().is_ok();

        Self {
            accession,
            name,
            aliases: Vec::new(),
            length,
            md5: None,
            assembly: GenomeAssembly::Custom,
            species,
            is_chromosome,
            chromosome: None,
        }
    }

    /// Get the chromosome identifier
    pub fn chromosome_id(&self) -> Option<&str> {
        if self.is_chromosome {
            Some(&self.name[3..])
        } else {
            self.chromosome.as_deref()
        }
    }

    /// Check if a position is valid for this reference
    pub fn is_valid_position(&self, pos: i32) -> bool {
        pos >= 0 && pos < self.length as i32
    }
}

/// Collection of reference sequences (e.g., all chromosomes in a genome)
#[derive(Debug, Clone)]
pub struct ReferenceCollection {
    /// Reference sequences
    references: Vec<ReferenceSequence>,
    /// Name to index mapping
    name_index: Vec<(String, usize)>,
}

impl ReferenceCollection {
    /// Create an empty collection
    pub fn new() -> Self {
        Self {
            references: Vec::new(),
            name_index: Vec::new(),
        }
    }

    /// Add a reference sequence
    pub fn add_reference(&mut self, ref_seq: ReferenceSequence) -> usize {
        let idx = self.references.len();
        self.references.push(ref_seq.clone());

        self.name_index.push((ref_seq.name.clone(), idx));

        for alias in &ref_seq.aliases {
            let alias_clone = alias.clone();
            self.name_index.push((alias_clone, idx));
        }

        if let Some(chrom) = ref_seq.chromosome_id() {
            self.name_index.push((chrom.to_string(), idx));
        }

        idx
    }

    /// Get reference by name
    pub fn get_by_name(&self, name: &str) -> Option<&ReferenceSequence> {
        self.name_index
            .iter()
            .find(|(n, _)| n == name)
            .map(|(_, idx)| &self.references[*idx])
    }

    /// Get reference by index
    pub fn get_by_index(&self, idx: usize) -> Option<&ReferenceSequence> {
        self.references.get(idx)
    }

    /// Get reference by accession
    pub fn get_by_accession(&self, accession: &str) -> Option<&ReferenceSequence> {
        self.references.iter().find(|r| r.accession == accession)
    }

    /// Get total number of references
    pub fn len(&self) -> usize {
        self.references.len()
    }

    /// Check if empty
    pub fn is_empty(&self) -> bool {
        self.references.is_empty()
    }

    /// Get all chromosome names
    pub fn chromosome_names(&self) -> Vec<&str> {
        self.references
            .iter()
            .filter(|r| r.is_chromosome)
            .map(|r| r.name.as_str())
            .collect()
    }

    /// Get total genome length
    pub fn total_length(&self) -> u64 {
        self.references.iter().map(|r| r.length as u64).sum()
    }
}

impl Default for ReferenceCollection {
    fn default() -> Self {
        Self::new()
    }
}

/// A genomic region specification
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GenomicRegion {
    /// Reference sequence name
    pub reference: String,
    /// Start position (0-based, inclusive)
    pub start: i32,
    /// End position (0-based, exclusive)
    pub end: i32,
}

impl GenomicRegion {
    /// Parse a region string (e.g., "chr1:100000-200000")
    pub fn parse(region: &str) -> Option<Self> {
        let parts: Vec<&str> = region.split(|c| c == ':' || c == '-').collect();
        if parts.len() != 3 {
            return None;
        }

        let reference = parts[0].to_string();
        let start = parts[1].parse().ok()?;
        let end = parts[2].parse().ok()?;

        if start >= end {
            return None;
        }

        Some(Self {
            reference,
            start,
            end,
        })
    }

    /// Create a region covering an entire reference
    pub fn full_reference(reference: String, length: i32) -> Self {
        Self {
            reference,
            start: 0,
            end: length,
        }
    }

    /// Get the length of the region
    pub fn len(&self) -> i32 {
        self.end - self.start
    }

    /// Check if a position is within this region
    pub fn contains(&self, pos: i32) -> bool {
        pos >= self.start && pos < self.end
    }

    /// Check if this region overlaps with another
    pub fn overlaps(&self, other: &GenomicRegion) -> bool {
        self.reference == other.reference && self.start < other.end && self.end > other.start
    }

    /// Expand the region by a given number of bases on each side
    pub fn expand(&self, left: i32, right: i32) -> Self {
        Self {
            reference: self.reference.clone(),
            start: (self.start - left).max(0),
            end: self.end + right,
        }
    }

    /// Get a string representation
    pub fn to_string(&self) -> String {
        format!("{}:{}-{}", self.reference, self.start, self.end)
    }
}

/// Strand direction
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strand {
    /// Forward strand (+)
    Forward,
    /// Reverse strand (-)
    Reverse,
}

impl Strand {
    /// Get the symbol for this strand
    pub fn symbol(&self) -> char {
        match self {
            Strand::Forward => '+',
            Strand::Reverse => '-',
        }
    }

    /// Parse from a character
    pub fn from_char(c: char) -> Option<Self> {
        match c {
            '+' | 'f' | 'F' => Some(Strand::Forward),
            '-' | 'r' | 'R' => Some(Strand::Reverse),
            _ => None,
        }
    }
}

/// A base position on a reference
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Position {
    /// Reference sequence name
    pub reference: String,
    /// Position (0-based)
    pub position: i32,
    /// Strand
    pub strand: Strand,
}

impl Position {
    /// Create a new position
    pub fn new(reference: String, position: i32, strand: Strand) -> Self {
        Self {
            reference,
            position,
            strand,
        }
    }

    /// Create a position on the forward strand
    pub fn forward(reference: String, position: i32) -> Self {
        Self::new(reference, position, Strand::Forward)
    }

    /// Create a position on the reverse strand
    pub fn reverse(reference: String, position: i32) -> Self {
        Self::new(reference, position, Strand::Reverse)
    }

    /// Get the 1-based position (common in genomics)
    pub fn one_based(&self) -> i32 {
        self.position + 1
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reference_sequence_chromosome_detection() {
        let chr1 = ReferenceSequence::new(
            "NC_000001.11".to_string(),
            "chr1".to_string(),
            248956422,
            "Homo sapiens".to_string(),
        );
        assert!(chr1.is_chromosome);

        let contig = ReferenceSequence::new(
            "NT_123456".to_string(),
            "KI270752.1".to_string(),
            100000,
            "Homo sapiens".to_string(),
        );
        assert!(!contig.is_chromosome);
    }

    #[test]
    fn test_reference_collection() {
        let mut collection = ReferenceCollection::new();

        collection.add_reference(ReferenceSequence::new(
            "NC_000001.11".to_string(),
            "chr1".to_string(),
            248956422,
            "Homo sapiens".to_string(),
        ));

        collection.add_reference(ReferenceSequence::new(
            "NC_000002.12".to_string(),
            "chr2".to_string(),
            242193529,
            "Homo sapiens".to_string(),
        ));

        assert_eq!(collection.len(), 2);

        assert!(collection.get_by_name("chr1").is_some());
        assert!(collection.get_by_name("NC_000001.11").is_some());
    }

    #[test]
    fn test_genomic_region_parsing() {
        let region = GenomicRegion::parse("chr1:100000-200000");
        assert!(region.is_some());

        let r = region.unwrap();
        assert_eq!(r.reference, "chr1");
        assert_eq!(r.start, 100000);
        assert_eq!(r.end, 200000);
        assert_eq!(r.len(), 100000);

        assert!(GenomicRegion::parse("invalid").is_none());
        assert!(GenomicRegion::parse("chr1:200000-100000").is_none());
    }

    #[test]
    fn test_region_expansion() {
        let region = GenomicRegion::parse("chr1:100-200").unwrap();
        let expanded = region.expand(10, 20);

        assert_eq!(expanded.reference, "chr1");
        assert_eq!(expanded.start, 90);
        assert_eq!(expanded.end, 220);
    }

    #[test]
    fn test_position_one_based() {
        let pos = Position::forward("chr1".to_string(), 0);
        assert_eq!(pos.one_based(), 1);

        let pos = Position::forward("chr1".to_string(), 99999);
        assert_eq!(pos.one_based(), 100000);
    }
}
