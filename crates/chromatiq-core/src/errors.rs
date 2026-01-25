//! Error types for Chromatiq Core

use thiserror::Error;

/// Core errors that can occur during genomic data processing
#[derive(Debug, Error)]
pub enum CoreError {
    /// File I/O errors
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),

    /// Invalid BAM file format
    #[error("Invalid BAM format: {0}")]
    InvalidFormat(&'static str),

    /// Invalid magic bytes in BAM header
    #[error("Invalid BAM magic bytes - file may be corrupted or not a BAM file")]
    InvalidMagic,

    /// Invalid compressed data
    #[error("Compression error: {0}")]
    Compression(String),

    /// Invalid reference sequence ID
    #[error("Invalid reference sequence ID: {0}")]
    InvalidReference(i32),

    /// Region out of bounds
    #[error("Region out of bounds: {0}")]
    OutOfBounds(String),

    /// Memory mapping errors
    #[error("Memory mapping error: {0}")]
    MmapError(String),

    /// Index file errors
    #[error("Index error: {0}")]
    IndexError(String),

    /// Unsupported feature
    #[error("Unsupported feature: {0}")]
    Unsupported(String),
}

/// Result type alias for Core operations
pub type CoreResult<T> = Result<T, CoreError>;
