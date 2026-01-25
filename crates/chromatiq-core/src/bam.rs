//! BAM file parsing and reading
//!
//! This module implements a high-performance BAM file parser.
//! BAM files are compressed binary formats containing aligned sequencing reads.

use alloc::string::String;
use alloc::vec::Vec;
use core::convert::TryFrom;
use memmap2::MmapOptions;

use crate::{CoreError, CoreResult};

/// BAM file magic bytes
const BAM_MAGIC: [u8; 4] = [b'B', b'A', b'M', 0x01];

/// Read alignment flags
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ReadFlags(pub u16);

impl ReadFlags {
    pub fn is_paired(self) -> bool {
        self.0 & 0x1 != 0
    }
    pub fn is_proper_pair(self) -> bool {
        self.0 & 0x2 != 0
    }
    pub fn is_unmapped(self) -> bool {
        self.0 & 0x4 != 0
    }
    pub fn is_mate_unmapped(self) -> bool {
        self.0 & 0x8 != 0
    }
    pub fn is_reverse(self) -> bool {
        self.0 & 0x10 != 0
    }
    pub fn is_mate_reverse(self) -> bool {
        self.0 & 0x20 != 0
    }
    pub fn is_first(self) -> bool {
        self.0 & 0x40 != 0
    }
    pub fn is_second(self) -> bool {
        self.0 & 0x80 != 0
    }
    pub fn is_secondary(self) -> bool {
        self.0 & 0x400 != 0
    }
    pub fn is_duplicate(self) -> bool {
        self.0 & 0x400 != 0
    }
    pub fn is_supplementary(self) -> bool {
        self.0 & 0x800 != 0
    }
}

/// CIGAR operation codes
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CigarOp {
    Match,
    SeqMatch,
    SeqMismatch,
    Deletion,
    Skipped,
    SoftClip,
    HardClip,
    Padding,
    Insertion,
}

impl TryFrom<u8> for CigarOp {
    type Error = CoreError;
    fn try_from(value: u8) -> Result<Self, Self::Error> {
        match value {
            b'M' => Ok(CigarOp::Match),
            b'=' => Ok(CigarOp::SeqMatch),
            b'X' => Ok(CigarOp::SeqMismatch),
            b'D' => Ok(CigarOp::Deletion),
            b'N' => Ok(CigarOp::Skipped),
            b'S' => Ok(CigarOp::SoftClip),
            b'H' => Ok(CigarOp::HardClip),
            b'P' => Ok(CigarOp::Padding),
            b'I' => Ok(CigarOp::Insertion),
            _ => Err(CoreError::InvalidFormat("Invalid CIGAR operation")),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CigarOperation {
    pub op: CigarOp,
    pub length: u32,
}

#[derive(Debug, Clone)]
pub struct BamRecord {
    pub ref_id: i32,
    pub pos: i32,
    pub name: String,
    pub mapq: u8,
    pub flags: ReadFlags,
    pub cigar: Vec<CigarOperation>,
    pub mate_ref_id: i32,
    pub mate_pos: i32,
    pub tlen: i32,
    pub sequence: Vec<u8>,
    pub quality: Vec<u8>,
    pub aux_data: Vec<AuxField>,
}

impl BamRecord {
    pub fn end_pos(&self) -> i32 {
        let mut end = self.pos;
        for cigar in &self.cigar {
            match cigar.op {
                CigarOp::Match
                | CigarOp::SeqMatch
                | CigarOp::SeqMismatch
                | CigarOp::Deletion
                | CigarOp::Skipped
                | CigarOp::Padding => {
                    end += cigar.length as i32;
                }
                _ => {}
            }
        }
        end
    }
    pub fn overlaps(&self, pos: i32) -> bool {
        pos >= self.pos && pos < self.end_pos()
    }
    pub fn seq_len(&self) -> usize {
        self.sequence.len()
    }
}

#[derive(Debug, Clone)]
pub struct AuxField {
    pub tag: [u8; 2],
    pub type_code: u8,
    pub value: AuxValue,
}

#[derive(Debug, Clone)]
pub enum AuxValue {
    Int8(i8),
    UInt8(u8),
    Int16(i16),
    UInt16(u16),
    Int32(i32),
    UInt32(u32),
    Float(f32),
    String(String),
    Hex(Vec<u8>),
}

#[derive(Debug, Clone)]
pub struct BamHeader {
    pub text: String,
    pub references: Vec<ReferenceInfo>,
    pub align_start: i32,
}

#[derive(Debug, Clone)]
pub struct ReferenceInfo {
    pub name: String,
    pub length: u32,
}

pub struct BamReader {
    mmap: memmap2::Mmap,
    pub header: BamHeader,
    block_start: usize,
    block_pos: usize,
}

impl BamReader {
    pub fn new(path: &std::path::Path) -> CoreResult<Self> {
        let file = std::fs::File::open(path).map_err(CoreError::Io)?;
        let mmap = unsafe {
            MmapOptions::new()
                .map(&file)
                .map_err(|e| CoreError::MmapError(e.to_string()))?
        };
        let mut reader = Self {
            mmap,
            header: BamHeader {
                text: String::new(),
                references: Vec::new(),
                align_start: 0,
            },
            block_start: 0,
            block_pos: 0,
        };
        reader.read_header()?;
        Ok(reader)
    }

    fn read_header(&mut self) -> CoreResult<()> {
        let data = &self.mmap;
        let mut pos = 0;
        if &data[pos..pos + 4] != BAM_MAGIC {
            return Err(CoreError::InvalidMagic);
        }
        pos += 4;
        let header_len = u32::from_le_bytes(data[pos..pos + 4].try_into().unwrap()) as usize;
        pos += 4;
        let header_text = core::str::from_utf8(&data[pos..pos + header_len])
            .map_err(|_| CoreError::InvalidFormat("Invalid header text encoding"))?
            .to_string();
        pos += header_len;
        let num_refs = i32::from_le_bytes(data[pos..pos + 4].try_into().unwrap()) as usize;
        pos += 4;
        let mut references = Vec::with_capacity(num_refs);
        for _ in 0..num_refs {
            let name_len = u32::from_le_bytes(data[pos..pos + 4].try_into().unwrap()) as usize;
            pos += 4;
            let name = core::str::from_utf8(&data[pos..pos + name_len - 1])
                .map_err(|_| CoreError::InvalidFormat("Invalid reference name encoding"))?
                .to_string();
            pos += name_len;
            let length = u32::from_le_bytes(data[pos..pos + 4].try_into().unwrap());
            pos += 4;
            references.push(ReferenceInfo { name, length });
        }
        self.header.text = header_text;
        self.header.references = references;
        self.header.align_start = pos as i32;
        self.block_start = pos;
        self.block_pos = pos;
        Ok(())
    }

    pub fn next_record(&mut self) -> Option<CoreResult<BamRecord>> {
        if self.block_pos >= self.mmap.len() {
            return None;
        }
        let data = &self.mmap[self.block_pos..];
        let block_size = u32::from_le_bytes(data[0..4].try_into().unwrap()) as usize;
        let record_data = &data[4..4 + block_size];
        let record = parse_bam_record(record_data);
        self.block_pos += 4 + block_size;
        Some(record)
    }

    pub fn get_reference(&self, ref_id: i32) -> Option<&ReferenceInfo> {
        if ref_id < 0 || ref_id >= self.header.references.len() as i32 {
            None
        } else {
            Some(&self.header.references[ref_id as usize])
        }
    }

    pub fn reference_names(&self) -> Vec<&str> {
        self.header
            .references
            .iter()
            .map(|r| r.name.as_str())
            .collect()
    }
}

fn parse_bam_record(input: &[u8]) -> CoreResult<BamRecord> {
    let mut pos = 0;
    let ref_id = i32::from_le_bytes(input[pos..pos + 4].try_into().unwrap());
    pos += 4;
    let pos_val = i32::from_le_bytes(input[pos..pos + 4].try_into().unwrap());
    pos += 4;
    let bin_mq_nl = u32::from_le_bytes(input[pos..pos + 4].try_into().unwrap());
    pos += 4;
    let bin = bin_mq_nl >> 16;
    let mapq = (bin_mq_nl >> 8 & 0xFF) as u8;
    let name_len = (bin_mq_nl & 0xFF) as usize;
    let flags_nc = u32::from_le_bytes(input[pos..pos + 4].try_into().unwrap());
    pos += 4;
    let flags = (flags_nc >> 16) as u16;
    let n_cigar_op = (flags_nc & 0xFFFF) as usize;
    let seq_len = i32::from_le_bytes(input[pos..pos + 4].try_into().unwrap()) as usize;
    pos += 4;
    let mate_ref_id = i32::from_le_bytes(input[pos..pos + 4].try_into().unwrap());
    pos += 4;
    let mate_pos = i32::from_le_bytes(input[pos..pos + 4].try_into().unwrap());
    pos += 4;
    let tlen = i32::from_le_bytes(input[pos..pos + 4].try_into().unwrap());
    pos += 4;
    let name_end = pos + name_len;
    let name = if name_len > 0 {
        core::str::from_utf8(&input[pos..name_end - 1])
            .map_err(|_| CoreError::InvalidFormat("Invalid read name encoding"))?
            .to_string()
    } else {
        String::new()
    };
    pos = name_end;
    let cigar_end = pos + n_cigar_op * 4;
    let cigar = parse_cigar(&input[pos..cigar_end], n_cigar_op)?;
    pos = cigar_end;
    let seq_end = pos + (seq_len + 1) / 2;
    let sequence = parse_sequence(&input[pos..seq_end], seq_len)?;
    pos = seq_end;
    let qual_end = pos + seq_len;
    let quality = input[pos..qual_end].to_vec();
    pos = qual_end;
    let aux_data = parse_aux_fields(&input[pos..])?;
    Ok(BamRecord {
        ref_id,
        pos: pos_val,
        name,
        mapq,
        flags: ReadFlags(flags),
        cigar,
        mate_ref_id,
        mate_pos,
        tlen,
        sequence,
        quality,
        aux_data,
    })
}

fn parse_cigar(input: &[u8], count: usize) -> CoreResult<Vec<CigarOperation>> {
    let mut cigar = Vec::with_capacity(count);
    let mut pos = 0;
    for _ in 0..count {
        if pos + 4 > input.len() {
            break;
        }
        let op_len = u32::from_le_bytes(input[pos..pos + 4].try_into().unwrap());
        pos += 4;
        let op_code = input[pos];
        pos += 1;
        let op = match op_code {
            b'M' => CigarOp::Match,
            b'I' => CigarOp::Insertion,
            b'D' => CigarOp::Deletion,
            b'N' => CigarOp::Skipped,
            b'S' => CigarOp::SoftClip,
            b'H' => CigarOp::HardClip,
            b'P' => CigarOp::Padding,
            b'=' => CigarOp::SeqMatch,
            b'X' => CigarOp::SeqMismatch,
            _ => CigarOp::Match,
        };
        cigar.push(CigarOperation { op, length: op_len });
    }
    Ok(cigar)
}

fn parse_sequence(input: &[u8], seq_len: usize) -> CoreResult<Vec<u8>> {
    let mut sequence = Vec::with_capacity(seq_len);
    for i in 0..seq_len {
        let byte = input[i / 2];
        sequence.push(if i % 2 == 0 { byte & 0x0F } else { byte >> 4 });
    }
    Ok(sequence)
}

fn parse_aux_fields(input: &[u8]) -> CoreResult<Vec<AuxField>> {
    let mut aux = Vec::new();
    let mut pos = 0;
    while pos + 2 < input.len() {
        let tag = [input[pos], input[pos + 1]];
        pos += 2;
        if pos >= input.len() {
            break;
        }
        let type_code = input[pos];
        pos += 1;
        let value = match type_code {
            b'Z' => {
                let end = pos;
                while end < input.len() && input[end] != 0 {
                    pos += 1;
                }
                AuxValue::String(
                    core::str::from_utf8(&input[end..pos])
                        .unwrap_or("")
                        .to_string(),
                )
            }
            _ => AuxValue::Int32(0),
        };
        aux.push(AuxField {
            tag,
            type_code,
            value,
        });
    }
    Ok(aux)
}

#[derive(Debug, Clone)]
pub struct ReferenceSequence {
    pub name: String,
    pub length: u32,
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_read_flags() {
        let flags = ReadFlags(0x43);
        assert!(flags.is_paired());
        assert!(flags.is_first());
        assert!(!flags.is_reverse());
    }
    #[test]
    fn test_cigar_ops() {
        assert_eq!(CigarOp::try_from(b'M').unwrap(), CigarOp::Match);
        assert!(CigarOp::try_from(b'?').is_err());
    }
}
