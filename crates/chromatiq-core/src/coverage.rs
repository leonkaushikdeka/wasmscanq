//! Coverage calculation algorithms
//!
//! This module provides efficient algorithms for calculating coverage depth
//! across genomic regions. Coverage represents how many times each base
//! has been sequenced in the dataset.

use alloc::vec::Vec;

use crate::{BamReader, BamRecord, CoreResult};

/// Coverage data for a genomic region
#[derive(Debug, Clone)]
pub struct Coverage {
    /// Reference sequence ID
    pub ref_id: i32,
    /// Reference name
    pub ref_name: String,
    /// Start position (0-based)
    pub start: i32,
    /// End position (exclusive)
    pub end: i32,
    /// Coverage values per position
    pub values: Vec<u32>,
    /// Minimum coverage depth
    pub min_depth: u32,
    /// Maximum coverage depth
    pub max_depth: u32,
    /// Mean coverage depth
    pub mean_depth: f64,
    /// Positions with zero coverage
    pub zero_coverage_positions: usize,
}

impl Coverage {
    /// Create a new empty coverage structure
    pub fn new(ref_id: i32, ref_name: String, start: i32, end: i32) -> Self {
        let len = (end - start) as usize;
        Self {
            ref_id,
            ref_name,
            start,
            end,
            values: vec![0; len],
            min_depth: 0,
            max_depth: 0,
            mean_depth: 0.0,
            zero_coverage_positions: 0,
        }
    }

    /// Get coverage at a specific position
    pub fn get(&self, pos: i32) -> Option<u32> {
        let idx = (pos - self.start) as usize;
        self.values.get(idx).copied()
    }

    /// Get coverage range for a region
    pub fn get_range(&self, start: i32, end: i32) -> Option<&[u32]> {
        let start_idx = (start - self.start) as usize;
        let end_idx = (end - self.start) as usize;
        self.values.get(start_idx..end_idx)
    }

    /// Calculate statistics after coverage is computed
    fn calculate_stats(&mut self) {
        if self.values.is_empty() {
            return;
        }

        let mut sum = 0u64;
        let mut zero_count = 0;

        for &val in &self.values {
            sum += val as u64;
            if val == 0 {
                zero_count += 1;
            }
        }

        self.min_depth = *self.values.iter().min().unwrap_or(&0);
        self.max_depth = *self.values.iter().max().unwrap_or(&0);
        self.mean_depth = sum as f64 / self.values.len() as f64;
        self.zero_coverage_positions = zero_count;
    }
}

/// Efficient coverage calculator using difference arrays
pub struct CoverageCalculator {
    /// Difference array for efficient range updates
    diff_array: Vec<i32>,
    /// Current region start
    start: i32,
    /// Current region end
    end: i32,
}

impl CoverageCalculator {
    /// Create a new coverage calculator for a region
    pub fn new(start: i32, end: i32) -> Self {
        let len = (end - start) as usize;
        Self {
            diff_array: vec![0; len + 1],
            start,
            end,
        }
    }

    /// Add coverage from a read
    pub fn add_read(&mut self, read: &BamRecord) {
        let read_start = read.pos;
        let read_end = read.end_pos();

        let overlap_start = core::cmp::max(read_start, self.start);
        let overlap_end = core::cmp::min(read_end, self.end);

        if overlap_start >= overlap_end {
            return;
        }

        let idx_start = (overlap_start - self.start) as usize;
        let idx_end = (overlap_end - self.start) as usize;

        self.diff_array[idx_start] += 1;
        if idx_end < self.diff_array.len() {
            self.diff_array[idx_end] -= 1;
        }
    }

    /// Add coverage from multiple reads in parallel
    pub fn add_reads_parallel(&mut self, reads: &[BamRecord]) {
        for read in reads {
            self.add_read(read);
        }
    }

    /// Convert difference array to actual coverage values
    pub fn finalize(self) -> Coverage {
        let values_len = self.diff_array.len() - 1;
        let mut values = Vec::with_capacity(values_len);

        let mut current = 0i32;
        for &diff in &self.diff_array[..values_len] {
            current += diff;
            values.push(current as u32);
        }

        let mut coverage = Coverage::new(0, String::new(), self.start, self.end);
        coverage.values = values;
        coverage.calculate_stats();

        coverage
    }
}

/// Calculate coverage for a region using a simple approach
pub fn calculate_coverage_simple(
    reader: &mut BamReader,
    ref_id: i32,
    start: i32,
    end: i32,
) -> CoreResult<Coverage> {
    let mut calculator = CoverageCalculator::new(start, end);

    while let Some(record_result) = reader.next_record() {
        let record = match record_result {
            Ok(r) => r,
            Err(_) => continue,
        };
        if record.ref_id != ref_id {
            continue;
        }
        calculator.add_read(&record);
    }

    let mut coverage = calculator.finalize();
    coverage.ref_id = ref_id;

    if let Some(ref_info) = reader.get_reference(ref_id) {
        coverage.ref_name = ref_info.name.clone();
    }

    Ok(coverage)
}

/// Calculate coverage with downsampling for large regions
pub fn calculate_coverage_downsampled(
    reader: &mut BamReader,
    ref_id: i32,
    start: i32,
    end: i32,
    max_points: usize,
) -> CoreResult<Coverage> {
    let region_size = end - start;
    if region_size <= max_points as i32 {
        return calculate_coverage_simple(reader, ref_id, start, end);
    }

    let bin_size = (region_size as usize + max_points - 1) / max_points;
    let mut binned_coverage = vec![0u64; max_points];

    while let Some(record_result) = reader.next_record() {
        let record = match record_result {
            Ok(r) => r,
            Err(_) => continue,
        };
        if record.ref_id != ref_id {
            continue;
        }

        let read_start = core::cmp::max(record.pos, start);
        let read_end = core::cmp::min(record.end_pos(), end);

        if read_start >= read_end {
            continue;
        }

        let start_bin = ((read_start - start) as usize) / bin_size;
        let end_bin = ((read_end - start) as usize).min(max_points - 1);

        for bin in start_bin..=end_bin {
            binned_coverage[bin] += 1;
        }
    }

    let mut coverage = Coverage::new(ref_id, String::new(), start, end);
    coverage.values = binned_coverage.iter().map(|&v| v as u32).collect();
    coverage.calculate_stats();

    if let Some(ref_info) = reader.get_reference(ref_id) {
        coverage.ref_name = ref_info.name.clone();
    }

    Ok(coverage)
}

/// High-resolution coverage calculation with base-level detail
pub struct HighResCoverage {
    pub ref_id: i32,
    pub ref_name: String,
    pub start: i32,
    pub end: i32,
    pub positions: Vec<i32>,
    pub depths: Vec<u32>,
}

impl HighResCoverage {
    /// Create a new high-resolution coverage tracker
    pub fn new(ref_id: i32, ref_name: String, start: i32, end: i32) -> Self {
        Self {
            ref_id,
            ref_name,
            start,
            end,
            positions: Vec::new(),
            depths: Vec::new(),
        }
    }

    /// Add a read to the coverage
    pub fn add_read(&mut self, read: &BamRecord) {
        let read_start = core::cmp::max(read.pos, self.start);
        let read_end = core::cmp::min(read.end_pos(), self.end);

        if read_start >= read_end {
            return;
        }

        for pos in read_start..read_end {
            if let Some(idx) = self.positions.iter().position(|&p| p == pos) {
                self.depths[idx] += 1;
            } else {
                self.positions.push(pos);
                self.depths.push(1);
            }
        }
    }

    /// Get coverage at a specific position using binary search
    pub fn get_coverage(&self, pos: i32) -> u32 {
        match self.positions.binary_search(&pos) {
            Ok(idx) => self.depths[idx],
            Err(_) => 0,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bam::{CigarOperation, ReadFlags};

    #[test]
    fn test_coverage_calculator_empty() {
        let calc = CoverageCalculator::new(0, 100);
        let coverage = calc.finalize();

        assert_eq!(coverage.values.len(), 100);
        assert!(coverage.values.iter().all(|&v| v == 0));
    }

    #[test]
    fn test_coverage_calculator_single_read() {
        let calc = CoverageCalculator::new(0, 100);

        let read = BamRecord {
            ref_id: 0,
            pos: 10,
            name: "read1".to_string(),
            mapq: 30,
            flags: ReadFlags(0),
            cigar: vec![CigarOperation {
                op: CigarOp::Match,
                length: 50,
            }],
            mate_ref_id: -1,
            mate_pos: 0,
            tlen: 0,
            sequence: vec![0; 50],
            quality: vec![30; 50],
            aux_data: Vec::new(),
        };

        calc.add_read(&read);
        let coverage = calc.finalize();

        for i in 10..60 {
            let idx = (i - 0) as usize;
            assert_eq!(coverage.values[idx], 1);
        }

        for i in [0, 1, 9, 60, 99] {
            let idx = (i - 0) as usize;
            assert_eq!(coverage.values[idx], 0);
        }
    }

    #[test]
    fn test_coverage_statistics() {
        let mut calc = CoverageCalculator::new(0, 100);

        for i in 0..10 {
            let read = BamRecord {
                ref_id: 0,
                pos: i * 10,
                name: format!("read{}", i),
                mapq: 30,
                flags: ReadFlags(0),
                cigar: vec![CigarOperation {
                    op: CigarOp::Match,
                    length: 10,
                }],
                mate_ref_id: -1,
                mate_pos: 0,
                tlen: 0,
                sequence: vec![0; 10],
                quality: vec![30; 10],
                aux_data: Vec::new(),
            };
            calc.add_read(&read);
        }

        let coverage = calc.finalize();

        assert!(coverage.max_depth > 0);
        assert!(coverage.mean_depth > 0.0);
    }
}
