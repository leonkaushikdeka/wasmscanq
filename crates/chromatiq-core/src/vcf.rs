//! VCF (Variant Call Format) file parsing and handling
//!
//! This module provides support for reading and parsing VCF files,
//! which contain variant calls from sequencing experiments.

use alloc::string::String;
use alloc::vec::Vec;

/// VCF record representing a single variant
#[derive(Debug, Clone)]
pub struct VcfRecord {
    /// Chromosome
    pub chrom: String,
    /// 1-based position
    pub pos: i32,
    /// Reference allele
    pub ref_allele: String,
    /// Alternate alleles
    pub alt_alleles: Vec<String>,
    /// Quality score
    pub qual: f32,
    /// Filter status (PASS, ., or list of failing filters)
    pub filter: String,
    /// INFO field values
    pub info: VcfInfo,
    /// Sample genotypes (if present)
    pub samples: Vec<VcfSample>,
}

/// VCF INFO field components
#[derive(Debug, Clone, Default)]
pub struct VcfInfo {
    /// Depth of coverage
    pub dp: Option<i32>,
    /// Allele depth for each allele
    pub ad: Option<Vec<i32>>,
    /// Genotype quality
    pub gq: Option<f32>,
    /// Genotype
    pub gt: Option<String>,
    /// Minor allele frequency
    pub maf: Option<f32>,
    /// Maximum likelihood estimate of variant frequency
    pub vaf: Option<f32>,
    /// Number of supporting reads
    pub support: Option<i32>,
    /// Strand bias
    pub strand_bias: Option<f32>,
}

/// Sample data from VCF
#[derive(Debug, Clone)]
pub struct VcfSample {
    pub name: String,
    pub genotype: String,
    pub gt_likelihoods: Vec<f32>,
    pub depth: i32,
    pub allele_depths: Vec<i32>,
    pub genotype_quality: f32,
}

/// VCF file header
#[derive(Debug, Clone)]
pub struct VcfHeader {
    pub fileformat: String,
    pub source: String,
    pub references: Vec<VcfReference>,
    pub contigs: Vec<VcfContig>,
    pub samples: Vec<String>,
    pub filters: Vec<VcfFilter>,
    pub infos: Vec<VcfInfoDef>,
    pub formats: Vec<VcfFormatDef>,
}

#[derive(Debug, Clone)]
pub struct VcfReference {
    pub id: String,
    pub url: String,
}

#[derive(Debug, Clone)]
pub struct VcfContig {
    pub id: String,
    pub length: i32,
}

#[derive(Debug, Clone)]
pub struct VcfFilter {
    pub id: String,
    pub description: String,
}

#[derive(Debug, Clone)]
pub struct VcfInfoDef {
    pub id: String,
    pub number: String,
    pub type_: String,
    pub description: String,
}

#[derive(Debug, Clone)]
pub struct VcfFormatDef {
    pub id: String,
    pub number: String,
    pub type_: String,
    pub description: String,
}

/// VCF parser
pub struct VcfReader {
    header: VcfHeader,
    records: Vec<VcfRecord>,
    current_pos: usize,
}

impl VcfReader {
    /// Create a new VCF reader from text content
    pub fn new(content: &str) -> Self {
        let mut lines = content.lines();
        let mut header = VcfHeader::default();
        let mut records = Vec::new();

        while let Some(line) = lines.next() {
            if line.starts_with("##") {
                parse_header_line(&line[2..], &mut header);
            } else if line.starts_with("#") {
                parse_column_header(&line[1..], &mut header);
            } else if !line.trim().is_empty() {
                if let Some(record) = parse_vcf_record(line, &header) {
                    records.push(record);
                }
            }
        }

        Self {
            header,
            records,
            current_pos: 0,
        }
    }

    /// Get file header
    pub fn header(&self) -> &VcfHeader {
        &self.header
    }

    /// Get next record
    pub fn next_record(&mut self) -> Option<&VcfRecord> {
        if self.current_pos < self.records.len() {
            let record = &self.records[self.current_pos];
            self.current_pos += 1;
            Some(record)
        } else {
            None
        }
    }

    /// Get all records
    pub fn records(&self) -> &[VcfRecord] {
        &self.records
    }

    /// Get records in a specific region
    pub fn records_in_region(&self, chrom: &str, start: i32, end: i32) -> Vec<&VcfRecord> {
        self.records
            .iter()
            .filter(|r| r.chrom == chrom && r.pos >= start && r.pos < end)
            .collect()
    }

    /// Get variant count
    pub fn variant_count(&self) -> usize {
        self.records.len()
    }

    /// Get samples
    pub fn samples(&self) -> &[String] {
        &self.header.samples
    }
}

impl Default for VcfHeader {
    fn default() -> Self {
        Self {
            fileformat: "VCFv4.2".to_string(),
            source: String::new(),
            references: Vec::new(),
            contigs: Vec::new(),
            samples: Vec::new(),
            filters: Vec::new(),
            infos: Vec::new(),
            formats: Vec::new(),
        }
    }
}

fn parse_header_line(line: &str, header: &mut VcfHeader) {
    if line.starts_with("FILEFORMAT=") {
        header.fileformat = line[11..].trim_matches('"').to_string();
    } else if line.starts_with("SOURCE=") {
        header.source = line[7..].trim_matches('"').to_string();
    } else if line.starts_with("REFERENCE=") {
        let parts: Vec<&str> = line[10..].splitn(2, '=').collect();
        if parts.len() == 2 {
            header.references.push(VcfReference {
                id: parts[0].to_string(),
                url: parts[1].trim_matches('"').to_string(),
            });
        }
    } else if line.starts_with("CONTIG=") {
        let contig_str = &line[7..];
        let mut id = String::new();
        let mut length = 0i32;

        for field in contig_str.split(',') {
            let parts: Vec<&str> = field.splitn(2, '=').collect();
            if parts.len() == 2 {
                match parts[0] {
                    "ID" => id = parts[1].trim_matches('"').to_string(),
                    "length" => length = parts[1].parse().unwrap_or(0),
                    _ => {}
                }
            }
        }
        if !id.is_empty() {
            header.contigs.push(VcfContig { id, length });
        }
    } else if line.starts_with("FILTER=") {
        let filter_str = &line[7..];
        let parts: Vec<&str> = filter_str.splitn(2, '=').collect();
        if parts.len() == 2 {
            header.filters.push(VcfFilter {
                id: parts[0].to_string(),
                description: parts[1].trim_matches('"').to_string(),
            });
        }
    } else if line.starts_with("INFO=") {
        let info_str = &line[5..];
        let mut id = String::new();
        let mut number = String::new();
        let mut type_ = String::new();
        let mut description = String::new();

        for field in info_str.split(',') {
            let parts: Vec<&str> = field.splitn(2, '=').collect();
            if parts.len() == 2 {
                match parts[0] {
                    "ID" => id = parts[1].to_string(),
                    "Number" => number = parts[1].to_string(),
                    "Type" => type_ = parts[1].to_string(),
                    "Description" => description = parts[1].trim_matches('"').to_string(),
                    _ => {}
                }
            }
        }
        header.infos.push(VcfInfoDef {
            id,
            number,
            type_,
            description,
        });
    }
}

fn parse_column_header(line: &str, header: &mut VcfHeader) {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() > 9 {
        for sample in &fields[9..] {
            header.samples.push(sample.to_string());
        }
    }
}

fn parse_vcf_record(line: &str, _header: &VcfHeader) -> Option<VcfRecord> {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 8 {
        return None;
    }

    let chrom = fields[0].to_string();
    let pos = fields[1].parse().ok()?;
    let id = fields[2].to_string();
    let ref_allele = fields[3].to_string();
    let alt_alleles: Vec<String> = if fields[4] == "." {
        Vec::new()
    } else {
        fields[4].split(',').map(|s| s.to_string()).collect()
    };
    let qual = fields[5].parse().unwrap_or(0.0);
    let filter = fields[6].to_string();

    let mut info = VcfInfo::default();
    if fields.len() > 7 && fields[7] != "." {
        parse_info_field(fields[7], &mut info);
    }

    let mut samples = Vec::new();
    if fields.len() > 9 {
        for (i, sample_data) in fields[9..].iter().enumerate() {
            if let Some(sample_name) = _header.samples.get(i) {
                if let Some(sample) = parse_sample(sample_data, sample_name) {
                    samples.push(sample);
                }
            }
        }
    }

    Some(VcfRecord {
        chrom,
        pos,
        ref_allele,
        alt_alleles,
        qual,
        filter,
        info,
        samples,
    })
}

fn parse_info_field(info_str: &str, info: &mut VcfInfo) {
    for field in info_str.split(';') {
        let parts: Vec<&str> = field.splitn(2, '=').collect();
        if parts.len() == 2 {
            match parts[0] {
                "DP" => info.dp = parts[1].parse().ok(),
                "AD" => {
                    info.ad = Some(parts[1].split(',').filter_map(|s| s.parse().ok()).collect());
                }
                "GQ" => info.gq = parts[1].parse().ok(),
                "GT" => info.gt = Some(parts[1].to_string()),
                "AF" => info.maf = parts[1].parse().ok(),
                "FREQ" | "VF" => info.vaf = parts[1].trim_end_matches('%').parse().ok(),
                "SUPPORT" | "SU" => info.support = parts[1].parse().ok(),
                "SB" => info.strand_bias = parts[1].parse().ok(),
                _ => {}
            }
        }
    }
}

fn parse_sample(sample_data: &str, _name: &str) -> Option<VcfSample> {
    let fields: Vec<&str> = sample_data.split(':').collect();
    if fields.is_empty() {
        return None;
    }

    let genotype = fields[0].to_string();
    let mut gt_likelihoods = Vec::new();
    let mut depth = 0i32;
    let mut allele_depths = Vec::new();
    let mut genotype_quality = 0.0f32;

    if fields.len() > 1 && fields[1] != "." {
        gt_likelihoods = fields[1]
            .split(',')
            .filter_map(|s| s.parse().ok())
            .collect();
    }

    if fields.len() > 2 && fields[2] != "." {
        if let Ok(d) = fields[2].parse() {
            depth = d;
        }
    }

    if fields.len() > 3 && fields[3] != "." {
        if let Ok(ad_str) = fields[3].parse::<String>() {
            allele_depths = ad_str.split(',').filter_map(|s| s.parse().ok()).collect();
        }
    }

    if fields.len() > 4 && fields[4] != "." {
        if let Ok(gq) = fields[4].parse() {
            genotype_quality = gq;
        }
    }

    Some(VcfSample {
        name: _name.to_string(),
        genotype,
        gt_likelihoods,
        depth,
        allele_depths,
        genotype_quality,
    })
}

/// Variant classification
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum VariantType {
    /// Single nucleotide polymorphism
    SNV,
    /// Multiple nucleotide polymorphism
    MNP,
    /// Insertion
    Insertion,
    /// Deletion
    Deletion,
    /// Structural variant
    SV,
    /// Unknown
    Unknown,
}

impl VcfRecord {
    /// Determine variant type
    pub fn variant_type(&self) -> VariantType {
        let ref_len = self.ref_allele.len();
        if self.alt_alleles.len() == 1 {
            let alt_len = self.alt_alleles[0].len();
            if ref_len == 1 && alt_len == 1 {
                VariantType::SNV
            } else if ref_len == alt_len {
                VariantType::MNP
            } else if alt_len > ref_len {
                VariantType::Insertion
            } else if alt_len < ref_len {
                VariantType::Deletion
            } else {
                VariantType::Unknown
            }
        } else {
            VariantType::Unknown
        }
    }

    /// Check if variant passed filters
    pub fn is_passed(&self) -> bool {
        self.filter == "PASS"
    }

    /// Get variant length
    pub fn variant_length(&self) -> i32 {
        let ref_len = self.ref_allele.len() as i32;
        let max_alt_len = self
            .alt_alleles
            .iter()
            .map(|s| s.len() as i32)
            .max()
            .unwrap_or(0);
        (max_alt_len - ref_len).abs()
    }

    /// Get allele frequency (if available)
    pub fn allele_frequency(&self) -> Option<f32> {
        self.info.vaf.or(self.info.maf)
    }

    /// Get supporting read count
    pub fn supporting_reads(&self) -> Option<i32> {
        self.info.support
    }
}

/// Summary statistics for VCF file
#[derive(Debug, Clone)]
pub struct VcfStats {
    pub total_variants: usize,
    pub snvs: usize,
    pub insertions: usize,
    pub deletions: usize,
    pub mnps: usize,
    pub passed_variants: usize,
    pub failed_variants: usize,
    pub high_quality_variants: usize,
    pub mean_quality: f32,
}

impl VcfStats {
    /// Calculate statistics from reader
    pub fn from_reader(reader: &VcfReader) -> Self {
        let mut snvs = 0;
        let mut insertions = 0;
        let mut deletions = 0;
        let mut mnps = 0;
        let mut passed = 0;
        let mut failed = 0;
        let mut high_quality = 0;
        let mut total_qual = 0.0f32;

        for record in reader.records() {
            match record.variant_type() {
                VariantType::SNV => snvs += 1,
                VariantType::Insertion => insertions += 1,
                VariantType::Deletion => deletions += 1,
                VariantType::MNP => mnps += 1,
                _ => {}
            }

            if record.is_passed() {
                passed += 1;
            } else {
                failed += 1;
            }

            if record.qual >= 30.0 {
                high_quality += 1;
            }

            total_qual += record.qual;
        }

        let total = reader.variant_count();
        Self {
            total_variants: total,
            snvs,
            insertions,
            deletions,
            mnps,
            passed_variants: passed,
            failed_variants: failed,
            high_quality_variants: high_quality,
            mean_quality: if total > 0 {
                total_qual / total as f32
            } else {
                0.0
            },
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_vcf_record_parsing() {
        let vcf_content = r#"##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##INFO=<ID=VT,Number=0,Type=String,Description="Variant Type">
##contig=<ID=chr1,length=249250621>
##contig=<ID=chr2,length=243199373>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	100000	.	A	G	50	PASS	DP=100;AF=0.25;VT=SNV
chr1	200000	.	AC	A	30	PASS	DP=50;AF=0.5;VT=DEL
chr1	300000	.	G	.	10	FAIL	DP=20
chr2	50000	.	A	T	60	PASS	DP=200;AF=0.75
"#;

        let reader = VcfReader::new(vcf_content);
        let records: Vec<_> = reader.records().iter().collect();

        assert_eq!(records.len(), 4);
        assert_eq!(records[0].chrom, "chr1");
        assert_eq!(records[0].pos, 100000);
        assert_eq!(records[0].ref_allele, "A");
        assert_eq!(records[0].alt_alleles[0], "G");
        assert!(records[0].is_passed());
        assert_eq!(records[0].info.dp, Some(100));
        assert_eq!(records[0].info.maf, Some(0.25));

        assert_eq!(records[1].variant_type(), VariantType::Deletion);
        assert_eq!(records[2].is_passed(), false);

        let stats = VcfStats::from_reader(&reader);
        assert_eq!(stats.total_variants, 4);
        assert_eq!(stats.snvs, 2);
        assert_eq!(stats.deletions, 1);
        assert_eq!(stats.passed_variants, 3);
        assert_eq!(stats.high_quality_variants, 3);
    }

    #[test]
    fn test_variant_type_detection() {
        let snv = VcfRecord {
            chrom: "chr1".to_string(),
            pos: 100,
            ref_allele: "A".to_string(),
            alt_alleles: vec!["G".to_string()],
            qual: 50.0,
            filter: "PASS".to_string(),
            info: VcfInfo::default(),
            samples: Vec::new(),
        };
        assert_eq!(snv.variant_type(), VariantType::SNV);

        let ins = VcfRecord {
            chrom: "chr1".to_string(),
            pos: 100,
            ref_allele: "A".to_string(),
            alt_alleles: vec!["ATGC".to_string()],
            qual: 50.0,
            filter: "PASS".to_string(),
            info: VcfInfo::default(),
            samples: Vec::new(),
        };
        assert_eq!(ins.variant_type(), VariantType::Insertion);
    }
}
