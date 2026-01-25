//! Chromatiq WASM Bindings
//!
//! This module provides WebAssembly bindings for the Chromatiq genomic
//! visualization engine, enabling high-performance genomic data processing
//! directly in web browsers.

#![cfg(target_arch = "wasm32")]

use serde::{Deserialize, Serialize};
use std::path::PathBuf;
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
    BamReader, BamRecord, CoreResult, Coverage, GenomicRegion, Pileup, PileupConfig,
    ReferenceCollection,
};

#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = BamFile)]
    pub type BamFile;
}

/// BAM file wrapper for WASM
#[wasm_bindgen]
pub struct BamFile {
    reader: Option<BamReader>,
    path: String,
    references: Vec<ReferenceInfo>,
    is_loaded: bool,
}

#[derive(Serialize, Deserialize)]
struct ReferenceInfo {
    name: String,
    length: u32,
}

#[wasm_bindgen]
impl BamFile {
    /// Create a new BAM file from a byte array
    #[wasm_bindgen(constructor)]
    pub fn new(data: &[u8]) -> Result<BamFile, JsValue> {
        console_log!("Loading BAM file from memory ({} bytes)", data.len());

        let temp_dir = std::env::temp_dir();
        let temp_path = temp_dir.join(format!("chromatiq_{}.bam", rand::random::<u64>()));

        std::fs::write(&temp_path, data)
            .map_err(|e| JsValue::from_str(&format!("Failed to write temp file: {}", e)))?;

        let reader = BamReader::new(&temp_path)
            .map_err(|e| JsValue::from_str(&format!("Failed to open BAM: {:?}", e)))?;

        let references = reader
            .reference_names()
            .iter()
            .map(|name| {
                let ref_info = reader.get_reference(0).cloned();
                ReferenceInfo {
                    name: name.to_string(),
                    length: reader.get_reference(0).map(|r| r.length).unwrap_or(0),
                }
            })
            .collect();

        Ok(Self {
            reader: Some(reader),
            path: temp_path.to_string_lossy().to_string(),
            references,
            is_loaded: true,
        })
    }

    /// Get list of reference sequences
    #[wasm_bindgen]
    pub fn get_references(&self) -> Result<JsValue, JsValue> {
        let refs: Vec<JsValue> = self
            .references
            .iter()
            .map(|r| serde_wasm_bindgen::to_value(r).unwrap())
            .collect();
        Ok(JsValue::from(&refs))
    }

    /// Get reference names as array
    #[wasm_bindgen]
    pub fn get_reference_names(&self) -> Vec<String> {
        self.references.iter().map(|r| r.name.clone()).collect()
    }

    /// Calculate coverage for a region
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

    /// Generate pileup for a region
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

    /// Get read count
    #[wasm_bindgen]
    pub fn get_read_count(&mut self) -> Result<u64, JsValue> {
        let reader = self
            .reader
            .as_mut()
            .ok_or_else(|| JsValue::from_str("BAM file not loaded"))?;

        let mut count = 0;
        while let Some(Some(_)) = reader.next_record().map(|r| r.ok()) {
            count += 1;
        }
        Ok(count)
    }
}

impl Drop for BamFile {
    fn drop(&mut self) {
        if self.is_loaded {
            let _ = std::fs::remove_file(&self.path);
        }
    }
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
}

#[wasm_bindgen]
impl ChromatiqEngine {
    /// Create a new Chromatiq engine
    #[wasm_bindgen(constructor)]
    pub fn new() -> Self {
        console_log!("Chromatiq Engine initialized");
        Self {
            bam_files: Vec::new(),
        }
    }

    /// Load a BAM file from byte array
    #[wasm_bindgen]
    pub fn load_bam(&mut self, data: &[u8]) -> Result<usize, JsValue> {
        let bam_file = BamFile::new(data)?;
        let idx = self.bam_files.len();
        self.bam_files.push(bam_file);
        Ok(idx)
    }

    /// Unload a BAM file
    #[wasm_bindgen]
    pub fn unload_bam(&mut self, index: usize) -> Result<(), JsValue> {
        if index >= self.bam_files.len() {
            return Err(JsValue::from_str("Invalid BAM file index"));
        }
        self.bam_files.remove(index);
        Ok(())
    }

    /// Get number of loaded BAM files
    #[wasm_bindgen]
    pub fn get_file_count(&self) -> usize {
        self.bam_files.len()
    }

    /// Calculate coverage for a region across all loaded files
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
}

#[derive(Serialize, Deserialize)]
struct CombinedCoverageData {
    ref_name: String,
    start: i32,
    end: i32,
    values: Vec<u32>,
}

/// Initialize panic hook for better error messages
#[wasm_bindgen(start)]
pub fn init() {
    console_error_panic_hook::set_once();
}

/// Get version information
#[wasm_bindgen]
pub fn get_version() -> String {
    env!("CARGO_PKG_VERSION").to_string()
}

/// Get engine information
#[wasm_bindgen]
pub fn get_info() -> String {
    format!(
        "Chromatiq {} - WASM Genomic Visualization Engine",
        env!("CARGO_PKG_VERSION")
    )
}
