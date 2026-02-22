/// 2-bit DNA encoding: A=0b00, C=0b01, G=0b10, T=0b11
///
/// Packs k-mers (k <= 32) into a single u64 for O(1) hashing and comparison.

use serde::{Deserialize, Serialize};

/// Encoding lookup: ASCII byte -> 2-bit value. Invalid bases map to 0xFF.
const ENCODE_TABLE: [u8; 256] = {
    let mut table = [0xFFu8; 256];
    table[b'A' as usize] = 0b00;
    table[b'a' as usize] = 0b00;
    table[b'C' as usize] = 0b01;
    table[b'c' as usize] = 0b01;
    table[b'G' as usize] = 0b10;
    table[b'g' as usize] = 0b10;
    table[b'T' as usize] = 0b11;
    table[b't' as usize] = 0b11;
    table
};

/// Decoding lookup: 2-bit value -> ASCII byte.
const DECODE_TABLE: [u8; 4] = [b'A', b'C', b'G', b'T'];

/// Complement lookup: 2-bit value -> complement 2-bit value.
/// A(00) <-> T(11), C(01) <-> G(10)
const COMPLEMENT_TABLE: [u8; 4] = [0b11, 0b10, 0b01, 0b00];

/// A k-mer packed into a u64 using 2-bit encoding.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug, Serialize, Deserialize)]
pub struct KmerU64 {
    pub encoded: u64,
    pub len: u8,
}

impl KmerU64 {
    /// Create a KmerU64 from a byte slice of ASCII DNA characters.
    /// Returns None if the slice contains non-ACGT characters or k > 32.
    #[inline]
    pub fn from_slice(seq: &[u8]) -> Option<Self> {
        let len = seq.len();
        if len == 0 || len > 32 {
            return None;
        }
        let mut encoded: u64 = 0;
        for &base in seq {
            let bits = ENCODE_TABLE[base as usize];
            if bits == 0xFF {
                return None;
            }
            encoded = (encoded << 2) | bits as u64;
        }
        Some(KmerU64 {
            encoded,
            len: len as u8,
        })
    }

    /// Decode back to a DNA string.
    pub fn to_string(&self) -> String {
        let mut result = Vec::with_capacity(self.len as usize);
        for i in (0..self.len).rev() {
            let bits = ((self.encoded >> (i * 2)) & 0b11) as usize;
            result.push(DECODE_TABLE[bits]);
        }
        String::from_utf8(result).unwrap()
    }

    /// Compute the reverse complement.
    pub fn reverse_complement(&self) -> KmerU64 {
        let mut rc: u64 = 0;
        let mut val = self.encoded;
        for _ in 0..self.len {
            let bits = (val & 0b11) as usize;
            rc = (rc << 2) | COMPLEMENT_TABLE[bits] as u64;
            val >>= 2;
        }
        KmerU64 {
            encoded: rc,
            len: self.len,
        }
    }

    /// Return the canonical form (lexicographically smaller of forward and reverse complement).
    pub fn canonical(&self) -> KmerU64 {
        let rc = self.reverse_complement();
        if self.encoded <= rc.encoded {
            *self
        } else {
            rc
        }
    }

    /// Split at position `pos` (0 < pos < len).
    /// Returns (left k-mer of length pos, right k-mer of length len-pos).
    #[inline]
    pub fn split_at(&self, pos: usize) -> (KmerU64, KmerU64) {
        debug_assert!(pos > 0 && pos < self.len as usize);
        let right_len = self.len as usize - pos;
        let right_mask = (1u64 << (right_len * 2)) - 1;
        let right_encoded = self.encoded & right_mask;
        let left_encoded = self.encoded >> (right_len * 2);
        (
            KmerU64 {
                encoded: left_encoded,
                len: pos as u8,
            },
            KmerU64 {
                encoded: right_encoded,
                len: right_len as u8,
            },
        )
    }
}

impl std::fmt::Display for KmerU64 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.to_string())
    }
}

/// Encode a single base character to its 2-bit representation.
/// Returns None for non-ACGT bases.
#[inline]
pub fn encode_base(base: u8) -> Option<u8> {
    let val = ENCODE_TABLE[base as usize];
    if val == 0xFF {
        None
    } else {
        Some(val)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_encode_decode_roundtrip() {
        let seqs = ["ACGT", "AAAA", "TTTT", "ACGTACGTACGTACGT", "A", "AT"];
        for s in seqs {
            let kmer = KmerU64::from_slice(s.as_bytes()).unwrap();
            assert_eq!(kmer.to_string(), s);
            assert_eq!(kmer.len as usize, s.len());
        }
    }

    #[test]
    fn test_reverse_complement() {
        let kmer = KmerU64::from_slice(b"ACGT").unwrap();
        let rc = kmer.reverse_complement();
        assert_eq!(rc.to_string(), "ACGT"); // ACGT is its own reverse complement

        let kmer2 = KmerU64::from_slice(b"AACG").unwrap();
        let rc2 = kmer2.reverse_complement();
        assert_eq!(rc2.to_string(), "CGTT");
    }

    #[test]
    fn test_canonical() {
        let kmer = KmerU64::from_slice(b"CGTT").unwrap();
        let canon = kmer.canonical();
        assert_eq!(canon.to_string(), "AACG");
    }

    #[test]
    fn test_split_at() {
        let kmer = KmerU64::from_slice(b"ACGT").unwrap();
        let (left, right) = kmer.split_at(2);
        assert_eq!(left.to_string(), "AC");
        assert_eq!(right.to_string(), "GT");

        let (left, right) = kmer.split_at(1);
        assert_eq!(left.to_string(), "A");
        assert_eq!(right.to_string(), "CGT");

        let (left, right) = kmer.split_at(3);
        assert_eq!(left.to_string(), "ACG");
        assert_eq!(right.to_string(), "T");
    }

    #[test]
    fn test_invalid_bases() {
        assert!(KmerU64::from_slice(b"ACGN").is_none());
        assert!(KmerU64::from_slice(b"").is_none());
    }

    #[test]
    fn test_case_insensitive() {
        let upper = KmerU64::from_slice(b"ACGT").unwrap();
        let lower = KmerU64::from_slice(b"acgt").unwrap();
        assert_eq!(upper.encoded, lower.encoded);
    }
}
