//! Indexed BAM file support with BAI index
//!
//! This module provides support for reading indexed BAM files
//! and efficiently fetching reads in a specific genomic region.

use alloc::string::String;
use alloc::vec::Vec;

use crate::bam::{BamReader, BamRecord, ReferenceInfo};
use crate::CoreResult;

/// BAI (BAM Index) header
#[derive(Debug, Clone)]
pub struct BaiIndex {
    /// File signature (BAI)
    pub signature: [u8; 4],
    /// Format version
    pub version: i32,
    /// Sort order (0=unknown, 1=coordinate, 2=queryname)
    pub sort_order: i32,
    /// Number of reference sequences
    pub n_refs: i32,
    /// Linear index for each reference
    pub refs: Vec<BaiReferenceIndex>,
}

/// Index for a single reference sequence
#[derive(Debug, Clone)]
pub struct BaiReferenceIndex {
    /// Reference name
    pub name: String,
    /// Reference length
    pub length: i64,
    /// Block count in auxiliary data
    pub n_blocks: i32,
    /// Block start positions
    pub block_starts: Vec<i64>,
    /// Block end positions
    pub block_ends: Vec<i64>,
    /// Number of chunks
    pub n_chunks: i32,
    /// Chunk start file offsets
    pub chunk_starts: Vec<i64>,
    /// Chunk end file offsets
    pub chunk_ends: Vec<i64>,
}

/// Indexed BAM reader that supports region queries
pub struct IndexedBamReader {
    bam_reader: BamReader,
    bai_index: Option<BaiIndex>,
    path: String,
}

impl IndexedBamReader {
    /// Open a BAM file with optional index
    pub fn new(bam_path: &str, bai_path: Option<&str>) -> CoreResult<Self> {
        let bam_reader = BamReader::new(&std::path::Path::new(bam_path))?;
        let bai_index = if let Some(bai) = bai_path {
            Some(BaiIndex::from_file(bai)?)
        } else {
            None
        };

        Ok(Self {
            bam_reader,
            bai_index,
            path: bam_path.to_string(),
        })
    }

    /// Get reference names
    pub fn reference_names(&self) -> Vec<&str> {
        self.bam_reader.reference_names()
    }

    /// Get reference info
    pub fn get_reference(&self, ref_id: i32) -> Option<&ReferenceInfo> {
        self.bam_reader.get_reference(ref_id)
    }

    /// Fetch reads in a region (requires BAI index for efficiency)
    pub fn fetch(&mut self, ref_name: &str, start: i32, end: i32) -> CoreResult<RegionIterator> {
        let ref_id = self
            .bam_reader
            .reference_names()
            .iter()
            .position(|n| n == &ref_name)
            .ok_or_else(|| crate::CoreError::InvalidReference(-1))? as i32;

        if let Some(ref index) = self.bai_index {
            if let Some(ref_idx) = index.refs.iter().find(|r| r.name == ref_name) {
                let chunk_offsets = self.find_chunks(ref_idx, start, end);

                if !chunk_offsets.is_empty() {
                    return Ok(RegionIterator {
                        reader: self,
                        ref_id,
                        ref_name: ref_name.to_string(),
                        start,
                        end,
                        chunk_offsets,
                        current_chunk: 0,
                        chunk_file_pos: 0,
                        last_alignment_start: -1,
                        record_buffer: Vec::new(),
                    });
                }
            }
        }

        Ok(RegionIterator {
            reader: self,
            ref_id,
            ref_name: ref_name.to_string(),
            start,
            end,
            chunk_offsets: Vec::new(),
            current_chunk: 0,
            chunk_file_pos: 0,
            last_alignment_start: -1,
            record_buffer: Vec::new(),
        })
    }

    fn find_chunks(&self, ref_idx: &BaiReferenceIndex, start: i32, end: i32) -> Vec<(i64, i64)> {
        let mut chunks = Vec::new();
        let bin_size: i64 = 16384;

        let start_bin = (start as i64 / bin_size) as usize;
        let end_bin = (end as i64 / bin_size) as usize;

        for (i, &chunk_start) in ref_idx.chunk_starts.iter().enumerate() {
            let chunk_end = ref_idx.chunk_ends[i];

            for bin_idx in start_bin..=end_bin.min(ref_idx.n_blocks as usize - 1) {
                if let Some((file_start, file_end)) = self.get_bin_data(ref_idx, bin_idx) {
                    if file_end > chunk_start && file_start < chunk_end {
                        chunks.push((chunk_start, chunk_end));
                        break;
                    }
                }
            }
        }

        chunks
    }

    fn get_bin_data(&self, _ref_idx: &BaiReferenceIndex, _bin_idx: usize) -> Option<(i64, i64)> {
        None
    }

    /// Get underlying BAM reader
    pub fn bam_reader(&mut self) -> &mut BamReader {
        &mut self.bam_reader
    }
}

impl BaiIndex {
    /// Load BAI index from file
    pub fn from_file(path: &str) -> CoreResult<Self> {
        let data = std::fs::read(path).map_err(|e| crate::CoreError::Io(e.into()))?;

        Self::from_bytes(&data)
    }

    /// Parse BAI index from bytes
    pub fn from_bytes(data: &[u8]) -> CoreResult<Self> {
        if data.len() < 16 {
            return Err(crate::CoreError::InvalidFormat("BAI file too small"));
        }

        let mut pos = 0;

        let signature = [data[0], data[1], data[2], data[3]];
        if signature != [b'B', b'A', b'I', 0x01] {
            return Err(crate::CoreError::InvalidMagic);
        }
        pos += 4;

        let version = i32::from_le_bytes(data[pos..pos + 4].try_into().unwrap());
        pos += 4;

        let sort_order = i32::from_le_bytes(data[pos..pos + 4].try_into().unwrap());
        pos += 4;

        let n_refs = i32::from_le_bytes(data[pos..pos + 4].try_into().unwrap());
        pos += 4;

        let mut refs = Vec::with_capacity(n_refs as usize);

        for _ in 0..n_refs {
            let name_len = i32::from_le_bytes(data[pos..pos + 4].try_into().unwrap());
            pos += 4;

            let name = if name_len > 0 {
                core::str::from_utf8(&data[pos..pos + name_len as usize - 1])
                    .map_err(|_| crate::CoreError::InvalidFormat("Invalid reference name"))?
                    .to_string()
            } else {
                String::new()
            };
            pos += name_len as usize;

            let length = i64::from_le_bytes(data[pos..pos + 8].try_into().unwrap());
            pos += 8;

            let n_blocks = i32::from_le_bytes(data[pos..pos + 4].try_into().unwrap());
            pos += 4;

            let mut block_starts = Vec::with_capacity(n_blocks as usize);
            let mut block_ends = Vec::with_capacity(n_blocks as usize);

            for _ in 0..n_blocks {
                let start = i64::from_le_bytes(data[pos..pos + 8].try_into().unwrap());
                pos += 8;
                let end = i64::from_le_bytes(data[pos..pos + 8].try_into().unwrap());
                pos += 8;
                block_starts.push(start);
                block_ends.push(end);
            }

            let n_chunks = i32::from_le_bytes(data[pos..pos + 4].try_into().unwrap());
            pos += 4;

            let mut chunk_starts = Vec::with_capacity(n_chunks as usize);
            let mut chunk_ends = Vec::with_capacity(n_chunks as usize);

            for _ in 0..n_chunks {
                let chunk_start = i64::from_le_bytes(data[pos..pos + 8].try_into().unwrap());
                pos += 8;
                let chunk_end = i64::from_le_bytes(data[pos..pos + 8].try_into().unwrap());
                pos += 8;
                chunk_starts.push(chunk_start);
                chunk_ends.push(chunk_end);
            }

            refs.push(BaiReferenceIndex {
                name,
                length,
                n_blocks,
                block_starts,
                block_ends,
                n_chunks,
                chunk_starts,
                chunk_ends,
            });
        }

        Ok(Self {
            signature,
            version,
            sort_order,
            n_refs,
            refs,
        })
    }
}

/// Iterator over reads in a genomic region
pub struct RegionIterator<'a> {
    reader: &'a mut IndexedBamReader,
    ref_id: i32,
    ref_name: String,
    start: i32,
    end: i32,
    chunk_offsets: Vec<(i64, i64)>,
    current_chunk: usize,
    chunk_file_pos: usize,
    last_alignment_start: i32,
    record_buffer: Vec<BamRecord>,
}

impl<'a> RegionIterator<'a> {
    /// Get the next read in the region
    pub fn next(&mut self) -> Option<CoreResult<BamRecord>> {
        if let Some(record) = self.record_buffer.pop() {
            return Some(Ok(record));
        }

        while let Some(record_result) = self.reader.bam_reader.next_record() {
            let record = match record_result {
                Ok(r) => r,
                Err(e) => return Some(Err(e)),
            };

            if record.ref_id != self.ref_id {
                continue;
            }

            if record.pos >= self.end {
                continue;
            }

            let record_end = record.end_pos();
            if record_end <= self.start {
                continue;
            }

            return Some(Ok(record));
        }

        None
    }

    /// Collect all reads into a vector
    pub fn collect(&mut self) -> CoreResult<Vec<BamRecord>> {
        let mut records = Vec::new();
        while let Some(record) = self.next() {
            records.push(record?);
        }
        Ok(records)
    }
}

/// Statistics for a region
#[derive(Debug, Clone)]
pub struct RegionStats {
    /// Total read count
    pub read_count: usize,
    /// Mapped read count
    pub mapped_count: usize,
    /// Unmapped read count
    pub unmapped_count: usize,
    /// Properly paired count
    pub properly_paired: usize,
    /// Duplicate count
    pub duplicate_count: usize,
    /// Mean mapping quality
    pub mean_mapq: f32,
    /// Total bases
    pub total_bases: usize,
    /// Mean read length
    pub mean_read_length: f32,
}

impl RegionStats {
    /// Calculate stats from records
    pub fn from_records(records: &[BamRecord]) -> Self {
        let mut read_count = 0;
        let mut mapped = 0;
        let mut properly_paired = 0;
        let mut duplicates = 0;
        let mut total_mapq = 0u64;
        let mut total_bases = 0usize;
        let mut total_length = 0u64;

        for record in records {
            read_count += 1;
            total_length += record.seq_len() as u64;
            total_bases += record.seq_len();

            if !record.flags.is_unmapped() {
                mapped += 1;
                total_mapq += record.mapq as u64;
            }

            if record.flags.is_proper_pair() {
                properly_paired += 1;
            }

            if record.flags.is_duplicate() {
                duplicates += 1;
            }
        }

        Self {
            read_count,
            mapped_count: mapped,
            unmapped_count: read_count - mapped,
            properly_paired,
            duplicate_count: duplicates,
            mean_mapq: if mapped > 0 {
                total_mapq as f32 / mapped as f32
            } else {
                0.0
            },
            total_bases,
            mean_read_length: if read_count > 0 {
                total_length as f32 / read_count as f32
            } else {
                0.0
            },
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bam::{CigarOp, CigarOperation, ReadFlags};

    #[test]
    fn test_bai_signature() {
        let valid_bai = [
            0x42, 0x41, 0x49, 0x01, // signature
            0x01, 0x00, 0x00, 0x00, // version
            0x01, 0x00, 0x00, 0x00, // sort order
            0x00, 0x00, 0x00, 0x00, // n_refs
        ];

        let index = BaiIndex::from_bytes(&valid_bai);
        assert!(index.is_ok());
        let idx = index.unwrap();
        assert_eq!(idx.signature, [b'B', b'A', b'I', 0x01]);
        assert_eq!(idx.version, 1);
        assert_eq!(idx.n_refs, 0);
    }

    #[test]
    fn test_region_stats() {
        let records = vec![
            BamRecord {
                ref_id: 0,
                pos: 100,
                name: "read1".to_string(),
                mapq: 30,
                flags: ReadFlags(99),
                cigar: vec![CigarOperation {
                    op: CigarOp::Match,
                    length: 50,
                }],
                mate_ref_id: 0,
                mate_pos: 200,
                tlen: 150,
                sequence: vec![b'A'; 50],
                quality: vec![30; 50],
                aux_data: Vec::new(),
            },
            BamRecord {
                ref_id: 0,
                pos: 150,
                name: "read2".to_string(),
                mapq: 20,
                flags: ReadFlags(73),
                cigar: vec![CigarOperation {
                    op: CigarOp::Match,
                    length: 75,
                }],
                mate_ref_id: -1,
                mate_pos: 0,
                tlen: 0,
                sequence: vec![b'G'; 75],
                quality: vec![25; 75],
                aux_data: Vec::new(),
            },
        ];

        let stats = RegionStats::from_records(&records);
        assert_eq!(stats.read_count, 2);
        assert_eq!(stats.mapped_count, 2);
        assert_eq!(stats.total_bases, 125);
        assert!((stats.mean_read_length - 62.5).abs() < 0.1);
    }
}
