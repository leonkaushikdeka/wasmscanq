//! Pileup generation algorithms
//!
//! Pileup is a columnar representation of aligned reads where each column
//! represents a genomic position with all overlapping reads displayed vertically.
//! This format is essential for variant calling and visualization.

use alloc::string::String;
use alloc::vec::Vec;

use crate::{bam::CigarOp, BamReader, BamRecord, CoreResult};

/// A single column in the pileup
#[derive(Debug, Clone)]
pub struct PileupColumn {
    /// Genomic position (0-based)
    pub position: i32,
    /// Reference base (if available)
    pub ref_base: Option<u8>,
    /// Reads overlapping this position
    pub reads: Vec<PileupRead>,
    /// Total depth at this position
    pub depth: usize,
    /// Consensus base (most common)
    pub consensus_base: Option<u8>,
    /// Consensus quality
    pub consensus_quality: f32,
}

/// A read in the pileup
#[derive(Debug, Clone)]
pub struct PileupRead {
    /// Read record
    pub record: BamRecord,
    /// Offset within the read (0-based)
    pub offset: usize,
    /// Base at this position
    pub base: u8,
    /// Quality score (Phred-scaled)
    pub quality: u8,
    /// Is this a deletion
    pub is_deletion: bool,
    /// Deletion length
    pub deletion_length: Option<u32>,
    /// Is this an insertion
    pub is_insertion: bool,
    /// Insertion sequence
    pub insertion_seq: Option<Vec<u8>>,
}

/// Complete pileup for a region
#[derive(Debug, Clone)]
pub struct Pileup {
    /// Reference sequence ID
    pub ref_id: i32,
    /// Reference name
    pub ref_name: String,
    /// Region start
    pub start: i32,
    /// Region end
    pub end: i32,
    /// Pileup columns
    pub columns: Vec<PileupColumn>,
    /// Total columns
    pub total_columns: usize,
    /// Maximum depth observed
    pub max_depth: usize,
}

/// Pileup configuration
#[derive(Debug, Clone)]
pub struct PileupConfig {
    /// Minimum mapping quality
    pub min_mapq: u8,
    /// Minimum base quality
    pub min_baseq: u8,
    /// Include duplicates
    pub include_duplicates: bool,
    /// Include secondary alignments
    pub include_secondary: bool,
    /// Include supplementary alignments
    pub include_supplementary: bool,
    /// Maximum reads per column
    pub max_depth: usize,
    /// Downsample probability (0.0 to 1.0)
    pub downsample: f32,
}

impl Default for PileupConfig {
    fn default() -> Self {
        Self {
            min_mapq: 0,
            min_baseq: 0,
            include_duplicates: false,
            include_secondary: false,
            include_supplementary: false,
            max_depth: 10000,
            downsample: 0.0,
        }
    }
}

/// Generator for pileup data
pub struct PileupGenerator {
    config: PileupConfig,
    current_reads: Vec<BamRecord>,
    current_pos: i32,
    end_pos: i32,
    ref_id: i32,
    sorted_reads: bool,
}

impl PileupGenerator {
    /// Create a new pileup generator
    pub fn new(ref_id: i32, start: i32, end: i32, config: PileupConfig) -> Self {
        Self {
            config,
            current_reads: Vec::new(),
            current_pos: start,
            end_pos: end,
            ref_id,
            sorted_reads: false,
        }
    }

    /// Set whether reads are coordinate-sorted
    pub fn set_sorted(&mut self, sorted: bool) {
        self.sorted_reads = sorted;
    }

    /// Add a read to the generator
    pub fn add_read(&mut self, read: BamRecord) {
        if read.ref_id != self.ref_id {
            return;
        }

        let read_end = read.end_pos();
        if read_end <= self.current_pos {
            return;
        }

        if self.current_reads.is_empty() {
            self.current_reads.push(read);
        } else {
            let insert_pos = self
                .current_reads
                .binary_search_by_key(&read.pos, |r| r.pos)
                .unwrap_or_else(|p| p);
            self.current_reads.insert(insert_pos, read);
        }
    }

    /// Generate the next column
    pub fn next_column(&mut self) -> Option<PileupColumn> {
        while self.current_pos < self.end_pos {
            let column = self.generate_column();
            if let Some(col) = &column {
                self.current_pos += 1;
                self.remove_finished_reads();
                return column;
            } else {
                self.current_pos += 1;
            }
        }
        None
    }

    fn generate_column(&mut self) -> Option<PileupColumn> {
        let mut reads_at_pos = Vec::new();

        for read in &self.current_reads {
            if !self.should_include_read(read) {
                continue;
            }

            let offset = (self.current_pos - read.pos) as usize;
            if offset >= read.sequence.len() {
                continue;
            }

            let (base, is_del, del_len, is_ins, ins_seq) = self.get_base_info(read, offset);

            reads_at_pos.push(PileupRead {
                record: read.clone(),
                offset,
                base,
                quality: read.quality.get(offset).copied().unwrap_or(0),
                is_deletion: is_del,
                deletion_length: del_len,
                is_insertion: is_ins,
                insertion_seq: ins_seq,
            });
        }

        if reads_at_pos.is_empty() {
            return None;
        }

        let max_depth = self.config.max_depth;
        let depth = reads_at_pos.len();
        if depth > max_depth {
            reads_at_pos.truncate(max_depth);
        }

        let (consensus_base, consensus_quality) = self.calculate_consensus(&reads_at_pos);

        Some(PileupColumn {
            position: self.current_pos,
            ref_base: None,
            reads: reads_at_pos,
            depth,
            consensus_base,
            consensus_quality,
        })
    }

    fn should_include_read(&self, read: &BamRecord) -> bool {
        if read.pos > self.current_pos {
            return false;
        }

        if read.end_pos() <= self.current_pos {
            return false;
        }

        if read.mapq < self.config.min_mapq {
            return false;
        }

        if !self.config.include_duplicates && read.flags.is_duplicate() {
            return false;
        }

        if !self.config.include_secondary && read.flags.is_secondary() {
            return false;
        }

        if !self.config.include_supplementary && read.flags.is_supplementary() {
            return false;
        }

        true
    }

    fn get_base_info(
        &self,
        read: &BamRecord,
        offset: usize,
    ) -> (u8, bool, Option<u32>, bool, Option<Vec<u8>>) {
        let mut base = b'N';
        let mut is_deletion = false;
        let mut deletion_length = None;
        let mut is_insertion = false;
        let mut insertion_seq = None;

        let mut cigar_offset = 0;
        let mut read_offset = 0;

        for cigar in &read.cigar {
            match cigar.op {
                CigarOp::Match | CigarOp::SeqMatch | CigarOp::SeqMismatch => {
                    let len = cigar.length as usize;
                    if offset >= cigar_offset && offset < cigar_offset + len {
                        if offset < read.sequence.len() {
                            base = read.sequence[offset];
                        }
                        break;
                    }
                    cigar_offset += len;
                }
                CigarOp::Deletion | CigarOp::Skipped => {
                    let len = cigar.length as usize;
                    if offset >= cigar_offset && offset < cigar_offset + len {
                        is_deletion = true;
                        deletion_length = Some(cigar.length);
                        break;
                    }
                    cigar_offset += len;
                }
                CigarOp::Insertion => {
                    let len = cigar.length as usize;
                    if offset >= cigar_offset && offset < cigar_offset + len {
                        is_insertion = true;
                        if offset < read.sequence.len() {
                            let ins_len = core::cmp::min(len, read.sequence.len() - offset);
                            insertion_seq = Some(read.sequence[offset..offset + ins_len].to_vec());
                        }
                        break;
                    }
                    cigar_offset += len;
                }
                _ => {}
            }
        }

        (
            base,
            is_deletion,
            deletion_length,
            is_insertion,
            insertion_seq,
        )
    }

    fn calculate_consensus(&self, reads: &[PileupRead]) -> (Option<u8>, f32) {
        if reads.is_empty() {
            return (None, 0.0);
        }

        let mut base_counts = [0u32; 256];
        let mut total_qual = 0u32;

        for read in reads {
            if read.is_deletion || read.quality < self.config.min_baseq {
                continue;
            }

            let base_idx = read.base as usize;
            base_counts[base_idx] += 1;
            total_qual += read.quality as u32;
        }

        let max_count = base_counts.iter().max().copied().unwrap_or(0);
        if max_count == 0 {
            return (None, 0.0);
        }

        let consensus_base = base_counts
            .iter()
            .enumerate()
            .find(|(_, &count)| count == max_count)
            .map(|(idx, _)| idx as u8)
            .unwrap_or(b'N');

        let consensus_quality = if total_qual > 0 {
            (total_qual as f32 / reads.len() as f32) - 33.0
        } else {
            0.0
        };

        (Some(consensus_base), consensus_quality)
    }

    fn remove_finished_reads(&mut self) {
        self.current_reads
            .retain(|read: &BamRecord| read.end_pos() > self.current_pos);
    }
}

/// Generate pileup for a region
pub fn generate_pileup(
    reader: &mut BamReader,
    ref_id: i32,
    start: i32,
    end: i32,
    config: Option<PileupConfig>,
) -> CoreResult<Pileup> {
    let config = config.unwrap_or_default();
    let mut generator = PileupGenerator::new(ref_id, start, end, config);

    while let Some(Some(record)) = reader.next_record().map(|r| r.ok()) {
        if record.ref_id == ref_id && record.pos < end && record.end_pos() > start {
            generator.add_read(record);
        }
    }

    let mut columns = Vec::new();
    let mut max_depth = 0;

    while let Some(column) = generator.next_column() {
        max_depth = core::cmp::max(max_depth, column.depth);
        columns.push(column);
    }

    let total_columns = columns.len();
    let mut pileup = Pileup {
        ref_id,
        ref_name: String::new(),
        start,
        end,
        columns,
        total_columns,
        max_depth,
    };

    if let Some(ref_info) = reader.get_reference(ref_id) {
        pileup.ref_name = ref_info.name.clone();
    }

    Ok(pileup)
}

/// Extract variant information from pileup
#[derive(Debug, Clone)]
pub struct VariantCandidate {
    pub position: i32,
    pub ref_base: u8,
    pub alt_base: u8,
    pub depth: usize,
    pub alt_depth: usize,
    pub frequency: f32,
    pub quality: f32,
}

impl VariantCandidate {
    /// Check if this is a high-confidence variant
    pub fn is_confident(&self, min_freq: f32, min_depth: usize, min_alt_count: usize) -> bool {
        self.alt_depth >= min_alt_count && self.depth >= min_depth && self.frequency >= min_freq
    }
}

/// Find variant candidates from pileup data
pub fn find_variants(pileup: &Pileup) -> Vec<VariantCandidate> {
    let mut variants = Vec::new();

    for column in &pileup.columns {
        if column.depth < 10 {
            continue;
        }

        let mut base_counts = [0u32; 256];
        let mut total_qual = 0u32;

        for read in &column.reads {
            if !read.is_deletion && read.quality >= 20 {
                base_counts[read.base as usize] += 1;
                total_qual += read.quality as u32;
            }
        }

        let ref_count = column
            .ref_base
            .map(|b| base_counts[b as usize])
            .unwrap_or(0);

        let total_variant_count = base_counts.iter().sum::<u32>() - ref_count;
        if total_variant_count < 3 {
            continue;
        }

        let max_alt_count = base_counts
            .iter()
            .enumerate()
            .filter(|(idx, _)| column.ref_base.map_or(true, |ref_b| *idx != ref_b as usize))
            .map(|(_, &count)| count)
            .max()
            .unwrap_or(0);

        if max_alt_count < 3 {
            continue;
        }

        let total = column.depth as f32;
        let alt_freq = max_alt_count as f32 / total;

        if alt_freq >= 0.05 && alt_freq <= 0.95 {
            let alt_base = base_counts
                .iter()
                .enumerate()
                .find(|(idx, &count)| {
                    column.ref_base.map_or(true, |ref_b| *idx != ref_b as usize)
                        && count == max_alt_count
                })
                .map(|(idx, _)| idx as u8)
                .unwrap_or(b'N');

            variants.push(VariantCandidate {
                position: column.position,
                ref_base: column.ref_base.unwrap_or(b'N'),
                alt_base,
                depth: column.depth,
                alt_depth: max_alt_count as usize,
                frequency: alt_freq,
                quality: total_qual as f32 / total,
            });
        }
    }

    variants
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bam::{CigarOperation, ReadFlags};

    fn create_test_read(pos: i32, cigar_len: u32) -> BamRecord {
        BamRecord {
            ref_id: 0,
            pos,
            name: format!("read_at_{}", pos),
            mapq: 30,
            flags: ReadFlags(0),
            cigar: vec![CigarOperation {
                op: CigarOp::Match,
                length: cigar_len,
            }],
            mate_ref_id: -1,
            mate_pos: 0,
            tlen: 0,
            sequence: vec![b'A'; cigar_len as usize],
            quality: vec![30; cigar_len as usize],
            aux_data: Vec::new(),
        }
    }

    #[test]
    fn test_pileup_generator_basic() {
        let config = PileupConfig::default();
        let mut gen = PileupGenerator::new(0, 0, 100, config);

        gen.add_read(create_test_read(10, 50));
        gen.add_read(create_test_read(20, 30));
        gen.add_read(create_test_read(15, 40));

        let column_10 = gen.next_column();
        assert!(column_10.is_some());
        assert_eq!(column_10.unwrap().position, 10);

        let column_15 =
            core::iter::from_fn(|| gen.next_column()).find(|column| column.position == 15);
        assert!(column_15.is_some());
        assert_eq!(column_15.unwrap().depth, 2);
    }

    #[test]
    fn test_pileup_generator_empty_region() {
        let config = PileupConfig::default();
        let mut gen = PileupGenerator::new(0, 0, 100, config);

        gen.add_read(create_test_read(200, 50));

        let column = gen.next_column();
        assert!(column.is_none());
    }

    #[test]
    fn test_pileup_config_defaults() {
        let config = PileupConfig::default();

        assert_eq!(config.min_mapq, 0);
        assert_eq!(config.min_baseq, 0);
        assert!(!config.include_duplicates);
        assert_eq!(config.max_depth, 10000);
    }

    #[test]
    fn test_variant_candidate_confidence() {
        let variant = VariantCandidate {
            position: 100,
            ref_base: b'A',
            alt_base: b'G',
            depth: 50,
            alt_depth: 25,
            frequency: 0.5,
            quality: 25.0,
        };

        assert!(variant.is_confident(0.1, 10, 3));
        assert!(!variant.is_confident(0.9, 10, 3));
        assert!(!variant.is_confident(0.1, 100, 3));
    }
}
