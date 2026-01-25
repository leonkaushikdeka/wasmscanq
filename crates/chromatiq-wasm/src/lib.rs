//! Chromatiq WASM Bindings
//!
//! This module provides WebAssembly bindings for the Chromatiq genomic
//! visualization engine, enabling high-performance genomic data processing
//! directly in web browsers.

#![cfg(target_arch = "wasm32")]

use serde::{Deserialize, Serialize};
use wasm_bindgen::prelude::*;

#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(js_namespace = console)]
    fn log(s: &str);
}

macro_rules! console_log {
    ($($t:tt)*) => (log(&format_args!($($t)*).to_string()))
}

use chromatiq_core::{
    bam::CigarOperation, bam::ReadFlags, BamReader, BamRecord, CoreResult, Coverage, GenomicRegion,
    Pileup, PileupConfig, ReferenceCollection, RegionStats, VariantType, VcfReader, VcfRecord,
    VcfStats,
};

#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = BamFile)]
    pub type BamFile;

    #[wasm_bindgen(typescript_type = VcfFile)]
    pub type VcfFile;
}

/// BAM file wrapper for WASM
#[wasm_bindgen]
pub struct BamFile {
    reader: Option<BamReader>,
    path: String,
    references: Vec<ReferenceInfo>,
    is_loaded: bool,
    read_count: u64,
}

#[derive(Serialize, Deserialize)]
struct ReferenceInfo {
    name: String,
    length: u32,
}

/// VCF file wrapper for WASM
#[wasm_bindgen]
pub struct VcfFile {
    reader: Option<VcfReader>,
    variants: Vec<VariantData>,
    stats: VcfStatsData,
}

#[derive(Serialize, Deserialize)]
struct VariantData {
    chrom: String,
    pos: i32,
    ref_allele: String,
    alt_alleles: Vec<String>,
    qual: f32,
    filter: String,
    variant_type: String,
    allele_frequency: Option<f32>,
    supporting_reads: Option<i32>,
    is_passed: bool,
}

#[derive(Serialize, Deserialize)]
struct VcfStatsData {
    total_variants: usize,
    snvs: usize,
    insertions: usize,
    deletions: usize,
    mnps: usize,
    passed_variants: usize,
    high_quality_variants: usize,
    mean_quality: f32,
}

#[wasm_bindgen]
impl BamFile {
    #[wasm_bindgen(constructor)]
    pub fn new(data: &[u8]) -> Result<BamFile, JsValue> {
        console_log!("Loading BAM file ({} bytes)", data.len());

        let temp_dir = std::env::temp_dir();
        let temp_path = temp_dir.join(format!("chromatiq_{}.bam", rand::random::<u64>()));

        std::fs::write(&temp_path, data)
            .map_err(|e| JsValue::from_str(&format!("Failed to write temp file: {}", e)))?;

        let reader = BamReader::new(&temp_path)
            .map_err(|e| JsValue::from_str(&format!("Failed to open BAM: {:?}", e)))?;

        let ref_names = reader.reference_names();
        let mut read_count = 0u64;

        for _ in reader.next_record() {
            read_count += 1;
        }

        let references: Vec<ReferenceInfo> = ref_names
            .iter()
            .enumerate()
            .map(|(i, name)| ReferenceInfo {
                name: name.to_string(),
                length: reader
                    .get_reference(i as i32)
                    .map(|r| r.length)
                    .unwrap_or(0),
            })
            .collect();

        Ok(Self {
            reader: Some(reader),
            path: temp_path.to_string_lossy().to_string(),
            references,
            is_loaded: true,
            read_count,
        })
    }

    #[wasm_bindgen]
    pub fn get_references(&self) -> Result<JsValue, JsValue> {
        let refs: Vec<JsValue> = self
            .references
            .iter()
            .map(|r| serde_wasm_bindgen::to_value(r).unwrap())
            .collect();
        Ok(JsValue::from(&refs))
    }

    #[wasm_bindgen]
    pub fn get_reference_names(&self) -> Vec<String> {
        self.references.iter().map(|r| r.name.clone()).collect()
    }

    #[wasm_bindgen]
    pub fn calculate_coverage(
        &mut self,
        ref_name: &str,
        start: i32,
        end: i32,
    ) -> Result<JsValue, JsValue> {
        let reader = self
            .reader
            .as_mut()
            .ok_or_else(|| JsValue::from_str("BAM file not loaded"))?;

        let ref_id = reader
            .reference_names()
            .iter()
            .position(|n| n == &ref_name)
            .ok_or_else(|| JsValue::from_str(&format!("Reference '{}' not found", ref_name)))?
            as i32;

        let coverage =
            chromatiq_core::coverage::calculate_coverage_simple(reader, ref_id, start, end)
                .map_err(|e| JsValue::from_str(&format!("Coverage error: {:?}", e)))?;

        let serialized = serde_wasm_bindgen::to_value(&CoverageData::from(coverage))
            .map_err(|e| JsValue::from_str(&format!("Serialization error: {:?}", e)))?;

        Ok(serialized)
    }

    #[wasm_bindgen]
    pub fn generate_pileup(
        &mut self,
        ref_name: &str,
        start: i32,
        end: i32,
        min_mapq: u8,
        min_baseq: u8,
    ) -> Result<JsValue, JsValue> {
        let reader = self
            .reader
            .as_mut()
            .ok_or_else(|| JsValue::from_str("BAM file not loaded"))?;

        let ref_id = reader
            .reference_names()
            .iter()
            .position(|n| n == &ref_name)
            .ok_or_else(|| JsValue::from_str(&format!("Reference '{}' not found", ref_name)))?
            as i32;

        let mut config = PileupConfig::default();
        config.min_mapq = min_mapq;
        config.min_baseq = min_baseq;

        let pileup =
            chromatiq_core::pileup::generate_pileup(reader, ref_id, start, end, Some(config))
                .map_err(|e| JsValue::from_str(&format!("Pileup error: {:?}", e)))?;

        let serialized = serde_wasm_bindgen::to_value(&PileupData::from(pileup))
            .map_err(|e| JsValue::from_str(&format!("Serialization error: {:?}", e)))?;

        Ok(serialized)
    }

    #[wasm_bindgen]
    pub fn get_read_count(&self) -> u64 {
        self.read_count
    }

    #[wasm_bindgen]
    pub fn get_region_stats(
        &mut self,
        ref_name: &str,
        start: i32,
        end: i32,
    ) -> Result<JsValue, JsValue> {
        let reader = self
            .reader
            .as_mut()
            .ok_or_else(|| JsValue::from_str("BAM file not loaded"))?;

        let ref_id = reader
            .reference_names()
            .iter()
            .position(|n| n == &ref_name)
            .ok_or_else(|| JsValue::from_str(&format!("Reference '{}' not found", ref_name)))?
            as i32;

        let mut records = Vec::new();
        while let Some(record_result) = reader.next_record() {
            let record = match record_result {
                Ok(r) => r,
                Err(_) => continue,
            };
            if record.ref_id != ref_id {
                continue;
            }
            if record.pos >= end {
                continue;
            }
            if record.end_pos() <= start {
                continue;
            }
            records.push(record);
        }

        let stats = RegionStats::from_records(&records);

        let stats_data = RegionStatsData {
            read_count: stats.read_count,
            mapped_count: stats.mapped_count,
            properly_paired: stats.properly_paired,
            duplicate_count: stats.duplicate_count,
            mean_mapq: stats.mean_mapq,
            mean_read_length: stats.mean_read_length,
        };

        serde_wasm_bindgen::to_value(&stats_data)
            .map_err(|e| JsValue::from_str(&format!("Serialization error: {:?}", e)))
    }
}

impl Drop for BamFile {
    fn drop(&mut self) {
        if self.is_loaded {
            let _ = std::fs::remove_file(&self.path);
        }
    }
}

#[wasm_bindgen]
impl VcfFile {
    #[wasm_bindgen(constructor)]
    pub fn new(data: &[u8]) -> Result<VcfFile, JsValue> {
        console_log!("Loading VCF file ({} bytes)", data.len());

        let content = String::from_utf8(data.to_vec())
            .map_err(|e| JsValue::from_str(&format!("Invalid UTF-8: {:?}", e)))?;

        let reader = VcfReader::new(&content);
        let stats = VcfStats::from_reader(&reader);

        let variants: Vec<VariantData> = reader
            .records()
            .iter()
            .map(|r| VariantData {
                chrom: r.chrom.clone(),
                pos: r.pos,
                ref_allele: r.ref_allele.clone(),
                alt_alleles: r.alt_alleles.clone(),
                qual: r.qual,
                filter: r.filter.clone(),
                variant_type: format!("{:?}", r.variant_type()),
                allele_frequency: r.allele_frequency(),
                supporting_reads: r.supporting_reads(),
                is_passed: r.is_passed(),
            })
            .collect();

        let stats_data = VcfStatsData {
            total_variants: stats.total_variants,
            snvs: stats.snvs,
            insertions: stats.insertions,
            deletions: stats.deletions,
            mnps: stats.mnps,
            passed_variants: stats.passed_variants,
            high_quality_variants: stats.high_quality_variants,
            mean_quality: stats.mean_quality,
        };

        Ok(Self {
            reader: Some(reader),
            variants,
            stats: stats_data,
        })
    }

    #[wasm_bindgen]
    pub fn get_variants(&self) -> Result<JsValue, JsValue> {
        serde_wasm_bindgen::to_value(&self.variants)
            .map_err(|e| JsValue::from_str(&format!("Serialization error: {:?}", e)))
    }

    #[wasm_bindgen]
    pub fn get_stats(&self) -> Result<JsValue, JsValue> {
        serde_wasm_bindgen::to_value(&self.stats)
            .map_err(|e| JsValue::from_str(&format!("Serialization error: {:?}", e)))
    }

    #[wasm_bindgen]
    pub fn get_variants_in_region(
        &self,
        chrom: &str,
        start: i32,
        end: i32,
    ) -> Result<JsValue, JsValue> {
        let filtered: Vec<&VariantData> = self
            .variants
            .iter()
            .filter(|v| v.chrom == chrom && v.pos >= start && v.pos < end)
            .collect();

        serde_wasm_bindgen::to_value(&filtered)
            .map_err(|e| JsValue::from_str(&format!("Serialization error: {:?}", e)))
    }

    #[wasm_bindgen]
    pub fn get_passed_variants(&self) -> Result<JsValue, JsValue> {
        let passed: Vec<&VariantData> = self.variants.iter().filter(|v| v.is_passed).collect();

        serde_wasm_bindgen::to_value(&passed)
            .map_err(|e| JsValue::from_str(&format!("Serialization error: {:?}", e)))
    }

    #[wasm_bindgen]
    pub fn get_high_quality_variants(&self, min_qual: f32) -> Result<JsValue, JsValue> {
        let hq: Vec<&VariantData> = self
            .variants
            .iter()
            .filter(|v| v.qual >= min_qual)
            .collect();

        serde_wasm_bindgen::to_value(&hq)
            .map_err(|e| JsValue::from_str(&format!("Serialization error: {:?}", e)))
    }
}

#[derive(Serialize, Deserialize)]
struct RegionStatsData {
    read_count: usize,
    mapped_count: usize,
    properly_paired: usize,
    duplicate_count: usize,
    mean_mapq: f32,
    mean_read_length: f32,
}

/// Coverage data for JavaScript
#[derive(Serialize, Deserialize)]
struct CoverageData {
    ref_id: i32,
    ref_name: String,
    start: i32,
    end: i32,
    values: Vec<u32>,
    min_depth: u32,
    max_depth: u32,
    mean_depth: f64,
    zero_coverage_positions: usize,
}

impl From<Coverage> for CoverageData {
    fn from(coverage: Coverage) -> Self {
        Self {
            ref_id: coverage.ref_id,
            ref_name: coverage.ref_name,
            start: coverage.start,
            end: coverage.end,
            values: coverage.values,
            min_depth: coverage.min_depth,
            max_depth: coverage.max_depth,
            mean_depth: coverage.mean_depth,
            zero_coverage_positions: coverage.zero_coverage_positions,
        }
    }
}

/// Pileup data for JavaScript
#[derive(Serialize, Deserialize)]
struct PileupData {
    ref_id: i32,
    ref_name: String,
    start: i32,
    end: i32,
    columns: Vec<PileupColumnData>,
    total_columns: usize,
    max_depth: usize,
}

#[derive(Serialize, Deserialize)]
struct PileupColumnData {
    position: i32,
    depth: usize,
    consensus_base: Option<u8>,
    consensus_quality: f32,
    reads: Vec<PileupReadData>,
}

#[derive(Serialize, Deserialize)]
struct PileupReadData {
    name: String,
    offset: usize,
    base: u8,
    quality: u8,
    is_deletion: bool,
    is_reverse: bool,
}

impl From<Pileup> for PileupData {
    fn from(pileup: Pileup) -> Self {
        Self {
            ref_id: pileup.ref_id,
            ref_name: pileup.ref_name,
            start: pileup.start,
            end: pileup.end,
            columns: pileup
                .columns
                .iter()
                .map(|c| PileupColumnData {
                    position: c.position,
                    depth: c.depth,
                    consensus_base: c.consensus_base,
                    consensus_quality: c.consensus_quality,
                    reads: c
                        .reads
                        .iter()
                        .map(|r| PileupReadData {
                            name: r.record.name.clone(),
                            offset: r.offset,
                            base: r.base,
                            quality: r.quality,
                            is_deletion: r.is_deletion,
                            is_reverse: r.record.flags.is_reverse(),
                        })
                        .collect(),
                })
                .collect(),
            total_columns: pileup.total_columns,
            max_depth: pileup.max_depth,
        }
    }
}

/// Chromatiq engine for genomic visualization
#[wasm_bindgen]
pub struct ChromatiqEngine {
    bam_files: Vec<BamFile>,
    vcf_files: Vec<VcfFile>,
}

#[wasm_bindgen]
impl ChromatiqEngine {
    #[wasm_bindgen(constructor)]
    pub fn new() -> Self {
        console_log!("Chromatiq Engine initialized");
        Self {
            bam_files: Vec::new(),
            vcf_files: Vec::new(),
        }
    }

    #[wasm_bindgen]
    pub fn load_bam(&mut self, data: &[u8]) -> Result<usize, JsValue> {
        let bam_file = BamFile::new(data)?;
        let idx = self.bam_files.len();
        self.bam_files.push(bam_file);
        Ok(idx)
    }

    #[wasm_bindgen]
    pub fn load_vcf(&mut self, data: &[u8]) -> Result<usize, JsValue> {
        let vcf_file = VcfFile::new(data)?;
        let idx = self.vcf_files.len();
        self.vcf_files.push(vcf_file);
        Ok(idx)
    }

    #[wasm_bindgen]
    pub fn unload_bam(&mut self, index: usize) -> Result<(), JsValue> {
        if index >= self.bam_files.len() {
            return Err(JsValue::from_str("Invalid BAM file index"));
        }
        self.bam_files.remove(index);
        Ok(())
    }

    #[wasm_bindgen]
    pub fn unload_vcf(&mut self, index: usize) -> Result<(), JsValue> {
        if index >= self.vcf_files.len() {
            return Err(JsValue::from_str("Invalid VCF file index"));
        }
        self.vcf_files.remove(index);
        Ok(())
    }

    #[wasm_bindgen]
    pub fn get_bam_file_count(&self) -> usize {
        self.bam_files.len()
    }

    #[wasm_bindgen]
    pub fn get_vcf_file_count(&self) -> usize {
        self.vcf_files.len()
    }

    #[wasm_bindgen]
    pub fn calculate_combined_coverage(
        &mut self,
        ref_name: &str,
        start: i32,
        end: i32,
    ) -> Result<JsValue, JsValue> {
        let mut combined_values = vec![0u32; (end - start) as usize];

        for bam in &mut self.bam_files {
            let coverage = bam.calculate_coverage(ref_name, start, end)?;

            let cov_data: CoverageData = serde_wasm_bindgen::from_value(coverage)
                .map_err(|e| JsValue::from_str(&format!("Deserialization error: {:?}", e)))?;

            for (i, &val) in cov_data.values.iter().enumerate() {
                if i < combined_values.len() {
                    combined_values[i] += val;
                }
            }
        }

        let combined = CombinedCoverageData {
            ref_name: ref_name.to_string(),
            start,
            end,
            values: combined_values,
        };

        serde_wasm_bindgen::to_value(&combined)
            .map_err(|e| JsValue::from_str(&format!("Serialization error: {:?}", e)))
    }

    #[wasm_bindgen]
    pub fn get_all_variants(&self) -> Result<JsValue, JsValue> {
        let all_variants: Vec<JsValue> = self
            .vcf_files
            .iter()
            .enumerate()
            .map(|(file_idx, vcf)| {
                let mut file_variants = vcf.variants.clone();
                for v in &mut file_variants {
                    v.chrom = format!("[{}] {}", file_idx, v.chrom);
                }
                serde_wasm_bindgen::to_value(&file_variants).unwrap_or(JsValue::NULL)
            })
            .collect();

        serde_wasm_bindgen::to_value(&all_variants)
            .map_err(|e| JsValue::from_str(&format!("Serialization error: {:?}", e)))
    }
}

#[derive(Serialize, Deserialize)]
struct CombinedCoverageData {
    ref_name: String,
    start: i32,
    end: i32,
    values: Vec<u32>,
}

#[wasm_bindgen(start)]
pub fn init() {
    console_error_panic_hook::set_once();
}

#[wasm_bindgen]
pub fn get_version() -> String {
    env!("CARGO_PKG_VERSION").to_string()
}

#[wasm_bindgen]
pub fn get_info() -> String {
    format!(
        "Chromatiq {} - WASM Genomic Visualization Engine with VCF support",
        env!("CARGO_PKG_VERSION")
    )
}

#[wasm_bindgen]
pub fn get_capabilities() -> String {
    serde_json::to_string(&Capabilities {
        bam_support: true,
        vcf_support: true,
        cram_support: false,
        indexed_bam: false,
        multi_track: true,
        offline_mode: true,
    })
    .unwrap_or_default()
}

#[derive(Serialize)]
struct Capabilities {
    bam_support: bool,
    vcf_support: bool,
    cram_support: bool,
    indexed_bam: bool,
    multi_track: bool,
    offline_mode: bool,
}
