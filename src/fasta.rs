/// Zero-copy memory-mapped FASTA reader.
///
/// Parses a FASTA file to find sequence regions, then provides direct
/// byte-slice access into the mmap. No sequence data is copied —
/// k-mer counting works directly on the mmap bytes, skipping newlines inline.

use anyhow::{Context, Result};
use memmap2::Mmap;
use std::fs::File;
use std::path::Path;

/// A region of sequence data within the memory-mapped file.
pub struct SequenceRegion {
    pub id: String,
    /// Byte offset of the first sequence character after the header line.
    pub start: usize,
    /// Byte offset past the last sequence character (up to next header or EOF).
    pub end: usize,
    /// Number of valid ACGT bases in this region.
    pub base_count: u64,
}

/// Zero-copy FASTA reader backed by a memory-mapped file.
pub struct FastaReader {
    mmap: Mmap,
    regions: Vec<SequenceRegion>,
    genome_length: u64,
}

impl FastaReader {
    /// Open and parse a FASTA file. Only parses headers and counts bases —
    /// sequence data stays in the mmap, never copied.
    pub fn from_path(path: &Path) -> Result<Self> {
        let file = File::open(path)
            .with_context(|| format!("Failed to open FASTA file: {}", path.display()))?;
        let mmap = unsafe { Mmap::map(&file) }
            .with_context(|| format!("Failed to mmap FASTA file: {}", path.display()))?;

        let (regions, genome_length) = parse_regions(&mmap);

        tracing::debug!(
            "Parsed {} sequences, {} bases from {}",
            regions.len(),
            genome_length,
            path.display()
        );

        Ok(FastaReader {
            mmap,
            regions,
            genome_length,
        })
    }

    /// All sequence regions.
    pub fn regions(&self) -> &[SequenceRegion] {
        &self.regions
    }

    /// Total valid bases across all sequences.
    pub fn genome_length(&self) -> u64 {
        self.genome_length
    }

    /// Raw bytes for a sequence region (includes newlines — caller must skip them).
    pub fn region_bytes(&self, region: &SequenceRegion) -> &[u8] {
        &self.mmap[region.start..region.end]
    }
}

/// Scan the mmap to find sequence regions and count valid bases.
fn parse_regions(data: &[u8]) -> (Vec<SequenceRegion>, u64) {
    let mut regions = Vec::new();
    let mut total_bases: u64 = 0;
    let mut i = 0;

    while i < data.len() {
        // Find next header line
        if data[i] != b'>' {
            i += 1;
            continue;
        }

        // Parse header: extract ID (up to first whitespace)
        let header_start = i + 1;
        while i < data.len() && data[i] != b'\n' {
            i += 1;
        }
        let header = &data[header_start..i];
        let id_end = header
            .iter()
            .position(|&b| b == b' ' || b == b'\t' || b == b'\r')
            .unwrap_or(header.len());
        let id = String::from_utf8_lossy(&header[..id_end]).to_string();

        // Skip past newline
        if i < data.len() {
            i += 1;
        }
        let seq_start = i;

        // Scan sequence data until next header or EOF, counting valid bases
        let mut base_count: u64 = 0;
        while i < data.len() && data[i] != b'>' {
            match data[i] {
                b'A' | b'a' | b'C' | b'c' | b'G' | b'g' | b'T' | b't' => base_count += 1,
                _ => {}
            }
            i += 1;
        }

        if base_count > 0 {
            total_bases += base_count;
            regions.push(SequenceRegion {
                id,
                start: seq_start,
                end: i,
                base_count,
            });
        }
    }

    (regions, total_bases)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_regions() {
        let data = b">seq1 description\nACGTACGT\nAAAA\n>seq2\nTTTTNNNNGGGG\n";
        let (regions, total) = parse_regions(data);

        assert_eq!(regions.len(), 2);
        assert_eq!(regions[0].id, "seq1");
        assert_eq!(regions[0].base_count, 12); // ACGTACGT + AAAA
        assert_eq!(regions[1].id, "seq2");
        assert_eq!(regions[1].base_count, 8); // TTTT + GGGG (N skipped)
        assert_eq!(total, 20);
    }

    #[test]
    fn test_parse_empty() {
        let (regions, total) = parse_regions(b"");
        assert!(regions.is_empty());
        assert_eq!(total, 0);
    }

    #[test]
    fn test_region_bytes_contain_newlines() {
        let data = b">seq1\nACGT\nAAAA\n";
        let (regions, _) = parse_regions(data);
        let region_data = &data[regions[0].start..regions[0].end];
        // Raw region includes the newline between lines
        assert_eq!(region_data, b"ACGT\nAAAA\n");
    }
}
