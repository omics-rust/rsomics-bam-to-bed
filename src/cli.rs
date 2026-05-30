use std::num::NonZero;
use std::path::PathBuf;

use clap::Parser;
use rsomics_common::{CommonFlags, Result, Tool, ToolMeta};
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};

use rsomics_bam_to_bed::{BamToBedOpts, bam_to_bed};

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(name = "rsomics-bam-to-bed", disable_help_flag = true)]
pub struct Cli {
    /// Input BAM file.
    #[arg(short = 'i', long = "input", required = true)]
    pub input: PathBuf,

    /// Split spliced alignments (N in CIGAR) into separate exon-block records.
    #[arg(long = "split")]
    pub split: bool,

    /// Use NM tag (edit distance) as BED score instead of mapping quality.
    #[arg(long = "ed")]
    pub use_edit_distance: bool,

    /// Append CIGAR string as a 7th column.
    #[arg(long = "cigar")]
    pub cigar: bool,

    #[command(flatten)]
    pub common: CommonFlags,
}

impl Tool for Cli {
    fn meta() -> ToolMeta {
        META
    }

    fn common(&self) -> &CommonFlags {
        &self.common
    }

    fn execute(self) -> Result<()> {
        let opts = BamToBedOpts {
            split: self.split,
            use_edit_distance: self.use_edit_distance,
            cigar: self.cigar,
        };
        let threads = self
            .common
            .threads
            .and_then(NonZero::new)
            .unwrap_or_else(|| {
                std::thread::available_parallelism().unwrap_or(NonZero::new(1).unwrap())
            });
        let stdout = std::io::stdout();
        let mut out = stdout.lock();
        bam_to_bed(&self.input, &mut out, &opts, threads)?;
        Ok(())
    }
}

pub const HELP: HelpSpec = HelpSpec {
    name: META.name,
    version: META.version,
    tagline: "Convert BAM alignments to BED6 format (bedtools bamtobed equivalent).",
    origin: Some(Origin {
        upstream: "bedtools",
        upstream_license: "MIT",
        our_license: "MIT OR Apache-2.0",
        paper_doi: Some("10.1093/bioinformatics/btq033"),
    }),
    usage_lines: &["-i <BAM> [OPTIONS]"],
    sections: &[Section {
        title: "OPTIONS",
        flags: &[
            FlagSpec {
                short: Some('i'),
                long: "input",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("Path"),
                required: true,
                default: None,
                description: "Input BAM file",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "split",
                aliases: &[],
                value: None,
                type_hint: Some("bool"),
                required: false,
                default: None,
                description: "Split spliced reads (N in CIGAR) into per-exon BED records",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "ed",
                aliases: &[],
                value: None,
                type_hint: Some("bool"),
                required: false,
                default: None,
                description: "Use edit distance (NM tag) as BED score instead of MAPQ",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "cigar",
                aliases: &[],
                value: None,
                type_hint: Some("bool"),
                required: false,
                default: None,
                description: "Append CIGAR string as 7th column",
                why_default: None,
            },
            FlagSpec {
                short: Some('h'),
                long: "help",
                aliases: &[],
                value: None,
                type_hint: Some("bool"),
                required: false,
                default: None,
                description: "Show this help",
                why_default: None,
            },
        ],
    }],
    examples: &[Example {
        description: "Convert a sorted BAM to BED",
        command: "rsomics-bam-to-bed -i alignments.bam",
    }],
    json_result_schema_doc: None,
};

#[cfg(test)]
mod tests {
    use clap::CommandFactory;

    #[test]
    fn cli_definition_is_valid() {
        super::Cli::command().debug_assert();
    }
}
