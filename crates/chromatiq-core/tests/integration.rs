use chromatiq_core::{CoverageCalculator, PileupConfig};

#[test]
fn public_api_smoke_test() {
    let coverage = CoverageCalculator::new(0, 10).finalize();

    assert_eq!(coverage.values.len(), 10);
    assert_eq!(PileupConfig::default().max_depth, 10000);
}
