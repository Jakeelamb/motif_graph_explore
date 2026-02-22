/// NCBI genome fetching via the Datasets API.
///
/// Downloads genome FASTA files by accession number (e.g. GCF_000005845.2).
/// Reads API key from .env file or NCBI_API_KEY environment variable.

use anyhow::{bail, Context, Result};
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::time::{Duration, Instant};

/// Configuration for NCBI API access.
pub struct NcbiConfig {
    pub api_key: Option<String>,
}

impl NcbiConfig {
    /// Load configuration from .env file or environment variable.
    pub fn load() -> Self {
        // Try .env file first
        let api_key = Self::read_env_file()
            .or_else(|| std::env::var("NCBI_API_KEY").ok());

        if api_key.is_some() {
            tracing::info!("NCBI API key loaded");
        } else {
            tracing::warn!("No NCBI API key found (requests will be rate-limited)");
        }

        NcbiConfig { api_key }
    }

    fn read_env_file() -> Option<String> {
        let content = std::fs::read_to_string(".env").ok()?;
        for line in content.lines() {
            let line = line.trim();
            if line.starts_with('#') || line.is_empty() {
                continue;
            }
            if let Some(rest) = line.strip_prefix("NCBI_API_KEY=") {
                let value = rest.trim_matches('"').trim_matches('\'');
                if !value.is_empty() {
                    return Some(value.to_string());
                }
            }
        }
        None
    }
}

/// Simple rate limiter for API requests.
pub struct RateLimiter {
    min_interval: Duration,
    last_request: Option<Instant>,
}

impl RateLimiter {
    /// Create a rate limiter.
    /// With API key: 10 requests/sec. Without: 3 requests/sec.
    pub fn new(has_api_key: bool) -> Self {
        let min_interval = if has_api_key {
            Duration::from_millis(100)
        } else {
            Duration::from_millis(334)
        };
        RateLimiter {
            min_interval,
            last_request: None,
        }
    }

    /// Wait if needed to respect rate limit.
    pub fn wait(&mut self) {
        if let Some(last) = self.last_request {
            let elapsed = last.elapsed();
            if elapsed < self.min_interval {
                std::thread::sleep(self.min_interval - elapsed);
            }
        }
        self.last_request = Some(Instant::now());
    }
}

/// Assembly summary information from NCBI.
#[derive(Debug)]
pub struct AssemblySummary {
    pub accession: String,
    pub organism: String,
    pub assembly_name: String,
    pub assembly_level: String,
}

/// Look up assembly information from NCBI.
pub fn lookup_assembly(accession: &str, config: &NcbiConfig) -> Result<AssemblySummary> {
    let url = format!(
        "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{}/dataset_report",
        accession
    );

    let mut request = ureq::get(&url)
        .set("Accept", "application/json");

    if let Some(ref key) = config.api_key {
        request = request.set("api-key", key);
    }

    let response = request.call()
        .with_context(|| format!("Failed to look up assembly: {}", accession))?;

    let body: serde_json::Value = serde_json::from_reader(response.into_reader())
        .context("Failed to parse NCBI response")?;

    // Parse the response
    let reports = body["reports"]
        .as_array()
        .context("No reports in NCBI response")?;

    if reports.is_empty() {
        bail!("No assembly found for accession: {}", accession);
    }

    let report = &reports[0];
    let organism = report["organism"]["organism_name"]
        .as_str()
        .unwrap_or("Unknown")
        .to_string();
    let assembly_name = report["assembly_info"]["assembly_name"]
        .as_str()
        .unwrap_or("Unknown")
        .to_string();
    let assembly_level = report["assembly_info"]["assembly_level"]
        .as_str()
        .unwrap_or("Unknown")
        .to_string();

    Ok(AssemblySummary {
        accession: accession.to_string(),
        organism,
        assembly_name,
        assembly_level,
    })
}

/// Download a genome FASTA file from NCBI by accession.
pub fn download_genome(accession: &str, output: &Path, config: &NcbiConfig) -> Result<PathBuf> {
    let url = format!(
        "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{}/download?include_annotation_type=GENOME_FASTA",
        accession
    );

    tracing::info!("Downloading genome: {}", accession);

    let mut request = ureq::get(&url);

    if let Some(ref key) = config.api_key {
        request = request.set("api-key", key);
    }

    let response = request.call()
        .with_context(|| format!("Failed to download genome: {}", accession))?;

    // Read response body into memory (ZIP file)
    let mut zip_data = Vec::new();
    response.into_reader().read_to_end(&mut zip_data)
        .context("Failed to read download response")?;

    // Extract FASTA from ZIP
    let cursor = std::io::Cursor::new(&zip_data);
    let mut archive = zip::ZipArchive::new(cursor)
        .context("Failed to open ZIP archive from NCBI response")?;

    // Find the .fna file in the archive
    let mut fna_name = None;
    for i in 0..archive.len() {
        let file = archive.by_index(i)?;
        let name = file.name().to_string();
        if name.ends_with(".fna") {
            fna_name = Some(name);
            break;
        }
    }

    let fna_name = fna_name.context("No .fna file found in NCBI ZIP archive")?;
    tracing::info!("Extracting: {}", fna_name);

    let mut fna_file = archive.by_name(&fna_name)
        .context("Failed to open .fna file in archive")?;

    // Ensure output directory exists
    if let Some(parent) = output.parent() {
        std::fs::create_dir_all(parent)?;
    }

    let mut out_file = std::fs::File::create(output)
        .with_context(|| format!("Failed to create output file: {}", output.display()))?;

    let mut buf = vec![0u8; 8192];
    loop {
        let n = fna_file.read(&mut buf)?;
        if n == 0 {
            break;
        }
        out_file.write_all(&buf[..n])?;
    }

    tracing::info!("Wrote FASTA to: {}", output.display());
    Ok(output.to_path_buf())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rate_limiter_spacing() {
        let mut limiter = RateLimiter::new(false);
        let start = Instant::now();
        limiter.wait(); // first call — no wait
        limiter.wait(); // second call — should wait ~334ms
        let elapsed = start.elapsed();
        assert!(elapsed >= Duration::from_millis(300), "Rate limiter should enforce minimum interval");
    }

    #[test]
    fn test_rate_limiter_with_api_key() {
        let mut limiter = RateLimiter::new(true);
        let start = Instant::now();
        limiter.wait();
        limiter.wait();
        let elapsed = start.elapsed();
        assert!(elapsed >= Duration::from_millis(80), "Rate limiter with API key should have shorter interval");
        assert!(elapsed < Duration::from_millis(300), "Rate limiter with API key should be faster than without");
    }

    #[test]
    fn test_config_from_env() {
        // Test that config loading doesn't panic
        std::env::set_var("NCBI_API_KEY", "test_key_12345");
        let config = NcbiConfig::load();
        assert!(config.api_key.is_some());
        std::env::remove_var("NCBI_API_KEY");
    }

    #[test]
    #[ignore] // Requires network access and valid API key
    fn test_download_e_coli() {
        let config = NcbiConfig::load();
        let tmp = std::env::temp_dir().join("test_ecoli.fna");
        let result = download_genome("GCF_000005845.2", &tmp, &config);
        assert!(result.is_ok());
        assert!(tmp.exists());
        let content = std::fs::read_to_string(&tmp).unwrap();
        assert!(content.starts_with('>'));
        std::fs::remove_file(tmp).ok();
    }
}
