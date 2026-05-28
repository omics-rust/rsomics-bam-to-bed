//! Convert BAM alignments to BED6 format.
//!
//! ## Algorithm
//!
//! For each mapped alignment, emit one BED6 record:
//! `chrom  start  end  name  score  strand`
//!
//! - `start` is 0-based (BED convention), sourced directly from the BAM `pos`
//!   field which is also 0-based.
//! - `end` = start + reference span of the alignment.
//! - `score` is the mapping quality (MAPQ) unless `use_edit_distance` is set,
//!   in which case the NM auxiliary tag value is used.
//! - `strand` is `+` for forward, `-` for reverse-complemented reads.
//!
//! With `split`, any spliced alignment (containing N CIGAR operations) is
//! broken into individual exon blocks — one BED record per block. The name and
//! score come from the read; blocks are non-overlapping and non-adjacent.
//!
//! Unmapped reads (FLAG & 0x4) are always skipped.
//!
//! ## Reference
//!
//! `BEDTools` bamtobed — Quinlan & Hall (2010). Bioinformatics 26(6): 841–842.
//! DOI: 10.1093/bioinformatics/btq033

use std::io::{BufWriter, Write};
use std::num::NonZero;
use std::path::Path;

use rsomics_bamio::raw::{self, RawRecord};
use rsomics_common::{Result, RsomicsError};

// BAM CIGAR op codes (packed low-nibble encoding).
const OP_MATCH: u8 = 0; // M
const OP_DEL: u8 = 2; // D
const OP_SKIP: u8 = 3; // N  (intron splice)
const OP_EQ: u8 = 7; // =
const OP_DIFF: u8 = 8; // X

/// Returns the number of reference bases consumed by `op`.
#[inline]
fn ref_len(op: u8, len: u32) -> u64 {
    match op {
        OP_MATCH | OP_DEL | OP_SKIP | OP_EQ | OP_DIFF => u64::from(len),
        _ => 0,
    }
}

#[derive(Default)]
pub struct BamToBedOpts {
    /// Split spliced alignments (N in CIGAR) into per-exon records.
    pub split: bool,
    /// Use NM tag value as BED score instead of MAPQ.
    pub use_edit_distance: bool,
    /// Append CIGAR string as a 7th column.
    pub cigar: bool,
}

/// A single BED record (reusable, reduces per-record allocation).
struct Bed6<'a> {
    chrom: &'a str,
    start: u64,
    end: u64,
    name: &'a [u8],
    score: u32,
    strand: u8,
    cigar: Option<&'a str>,
}

impl Bed6<'_> {
    fn write(&self, out: &mut impl Write) -> Result<()> {
        let mut ib = itoa::Buffer::new();
        out.write_all(self.chrom.as_bytes())
            .map_err(RsomicsError::Io)?;
        out.write_all(b"\t").map_err(RsomicsError::Io)?;
        out.write_all(ib.format(self.start).as_bytes())
            .map_err(RsomicsError::Io)?;
        out.write_all(b"\t").map_err(RsomicsError::Io)?;
        out.write_all(ib.format(self.end).as_bytes())
            .map_err(RsomicsError::Io)?;
        out.write_all(b"\t").map_err(RsomicsError::Io)?;
        out.write_all(self.name).map_err(RsomicsError::Io)?;
        out.write_all(b"\t").map_err(RsomicsError::Io)?;
        out.write_all(ib.format(self.score).as_bytes())
            .map_err(RsomicsError::Io)?;
        out.write_all(b"\t").map_err(RsomicsError::Io)?;
        out.write_all(&[self.strand]).map_err(RsomicsError::Io)?;
        if let Some(c) = self.cigar {
            out.write_all(b"\t").map_err(RsomicsError::Io)?;
            out.write_all(c.as_bytes()).map_err(RsomicsError::Io)?;
        }
        out.write_all(b"\n").map_err(RsomicsError::Io)?;
        Ok(())
    }
}

pub fn bam_to_bed(
    input: &Path,
    output: &mut dyn Write,
    opts: &BamToBedOpts,
    workers: NonZero<usize>,
) -> Result<u64> {
    let mut reader = rsomics_bamio::open_with_workers(input, workers)?;
    let header = reader.read_header().map_err(RsomicsError::Io)?;
    let ref_names: Vec<String> = header
        .reference_sequences()
        .keys()
        .map(ToString::to_string)
        .collect();

    let mut out = BufWriter::with_capacity(256 * 1024, output);
    let mut record = RawRecord::default();
    let mut count: u64 = 0;

    while raw::read_record(reader.get_mut(), &mut record)? != 0 {
        let flags = record.flags();

        // Skip unmapped reads.
        if flags & 0x4 != 0 {
            continue;
        }

        let tid = record.reference_sequence_id();
        let pos0 = record.alignment_start(); // 0-based BAM pos
        if tid < 0 || pos0 < 0 {
            continue;
        }

        let chrom = usize::try_from(tid)
            .ok()
            .and_then(|i| ref_names.get(i))
            .map_or("*", String::as_str);
        let name = build_name(record.name(), flags);
        let strand = if flags & 0x10 != 0 { b'-' } else { b'+' };
        let score: u32 = if opts.use_edit_distance {
            nm_tag(&record)
        } else {
            u32::from(record.mapping_quality())
        };
        let cigar_str: Option<String> = if opts.cigar {
            Some(cigar_to_string(record.cigar_ops()))
        } else {
            None
        };

        if opts.split {
            let ctx = SplitContext {
                chrom,
                name: &name,
                score,
                strand,
                cigar_str: cigar_str.as_deref(),
            };
            count += emit_split_blocks(&record, &ctx, u64::try_from(pos0).unwrap_or(0), &mut out)?;
        } else {
            let ref_span: u64 = record.cigar_ops().map(|(op, len)| ref_len(op, len)).sum();
            let start = u64::try_from(pos0).unwrap_or(0);
            let end = start + ref_span;
            Bed6 {
                chrom,
                start,
                end,
                name: &name,
                score,
                strand,
                cigar: cigar_str.as_deref(),
            }
            .write(&mut out)?;
            count += 1;
        }
    }

    out.flush().map_err(RsomicsError::Io)?;
    Ok(count)
}

struct SplitContext<'a> {
    chrom: &'a str,
    name: &'a [u8],
    score: u32,
    strand: u8,
    cigar_str: Option<&'a str>,
}

/// Emit one BED record per exon block (CIGAR runs split at N ops).
fn emit_split_blocks(
    record: &RawRecord,
    ctx: &SplitContext<'_>,
    pos0: u64,
    out: &mut impl Write,
) -> Result<u64> {
    let SplitContext {
        chrom,
        name,
        score,
        strand,
        cigar_str,
    } = ctx;
    let mut ref_cursor = pos0;
    let mut block_start = pos0;
    let mut in_block = false;
    let mut count: u64 = 0;

    for (op, len) in record.cigar_ops() {
        let rlen = u64::from(len);
        match op {
            OP_SKIP => {
                // N: splice junction — emit the accumulated exon block.
                if in_block {
                    Bed6 {
                        chrom,
                        start: block_start,
                        end: ref_cursor,
                        name,
                        score: *score,
                        strand: *strand,
                        cigar: None,
                    }
                    .write(out)?;
                    count += 1;
                    in_block = false;
                }
                ref_cursor += rlen;
                block_start = ref_cursor;
            }
            OP_MATCH | OP_DEL | OP_EQ | OP_DIFF => {
                if !in_block {
                    block_start = ref_cursor;
                    in_block = true;
                }
                ref_cursor += rlen;
            }
            _ => {}
        }
    }

    // Emit the final block.
    if in_block {
        Bed6 {
            chrom,
            start: block_start,
            end: ref_cursor,
            name,
            score: *score,
            strand: *strand,
            cigar: *cigar_str,
        }
        .write(out)?;
        count += 1;
    }

    Ok(count)
}

/// Append `/1` or `/2` mate suffix to paired reads, matching bedtools bamtobed.
///
/// bedtools bamtobed emits the query name with `/1` for the first segment and
/// `/2` for the second segment of a paired-end read. Single-end reads and reads
/// already carrying a suffix (rare) get no modification.
fn build_name(raw: &[u8], flags: u16) -> Vec<u8> {
    let is_paired = flags & 0x1 != 0;
    if !is_paired {
        return raw.to_vec();
    }
    let suffix: &[u8] = if flags & 0x40 != 0 { b"/1" } else { b"/2" };
    let mut name = Vec::with_capacity(raw.len() + 2);
    name.extend_from_slice(raw);
    name.extend_from_slice(suffix);
    name
}

/// Read the NM (edit distance) auxiliary tag, returning 0 if absent.
fn nm_tag(record: &RawRecord) -> u32 {
    record
        .aux_value(*b"NM")
        .and_then(|val| {
            let type_code = record.aux_type(*b"NM")?;
            parse_int_aux(type_code, val)
        })
        .unwrap_or(0)
}

fn parse_int_aux(type_code: u8, val: &[u8]) -> Option<u32> {
    match type_code {
        // Signed byte: reinterpret bits then cast to u32 (edit distance is always >= 0).
        b'c' => val
            .first()
            .map(|&v| i32::from(v.cast_signed()).cast_unsigned()),
        b'C' => val.first().copied().map(u32::from),
        b's' => val
            .get(..2)
            .map(|b| i32::from(i16::from_le_bytes([b[0], b[1]])).cast_unsigned()),
        b'S' => val
            .get(..2)
            .map(|b| u32::from(u16::from_le_bytes([b[0], b[1]]))),
        b'i' => val
            .get(..4)
            .map(|b| i32::from_le_bytes([b[0], b[1], b[2], b[3]]).cast_unsigned()),
        b'I' => val
            .get(..4)
            .map(|b| u32::from_le_bytes([b[0], b[1], b[2], b[3]])),
        _ => None,
    }
}

const CIGAR_CHARS: &[u8] = b"MIDNSHP=X";

fn cigar_to_string(ops: impl Iterator<Item = (u8, u32)>) -> String {
    let mut s = String::new();
    let mut ib = itoa::Buffer::new();
    for (op, len) in ops {
        s.push_str(ib.format(len));
        s.push(CIGAR_CHARS.get(op as usize).copied().unwrap_or(b'?') as char);
    }
    s
}
