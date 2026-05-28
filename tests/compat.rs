/// Compatibility tests: rsomics-bam-to-bed output must match bedtools bamtobed.
use std::path::PathBuf;
use std::process::Command;

fn bin_path() -> &'static str {
    env!("CARGO_BIN_EXE_rsomics-bam-to-bed")
}

fn golden(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

fn bedtools_available() -> bool {
    Command::new("bedtools")
        .arg("--version")
        .output()
        .is_ok_and(|o| o.status.success())
}

fn run_ours(bam: &str) -> String {
    let out = Command::new(bin_path())
        .args(["-i", bam])
        .output()
        .expect("failed to run rsomics-bam-to-bed");
    assert!(
        out.status.success(),
        "exit {:?}\nstderr: {}",
        out.status,
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8(out.stdout).unwrap()
}

fn run_bedtools(bam: &str) -> String {
    let out = Command::new("bedtools")
        .args(["bamtobed", "-i", bam])
        .output()
        .expect("failed to run bedtools bamtobed");
    assert!(
        out.status.success(),
        "bedtools exit {:?}\nstderr: {}",
        out.status,
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8(out.stdout).unwrap()
}

#[test]
fn golden_output_matches() {
    let bam = golden("small.bam");
    let got = run_ours(bam.to_str().unwrap());
    let expected = std::fs::read_to_string(golden("expected.bed")).unwrap();
    assert_eq!(
        got, expected,
        "output diverged from golden\ngot:\n{got}\nexpected:\n{expected}"
    );
}

#[test]
fn output_matches_bedtools() {
    if !bedtools_available() {
        eprintln!("skip: bedtools not found");
        return;
    }
    let bam = golden("small.bam");
    let ours = run_ours(bam.to_str().unwrap());
    let bt = run_bedtools(bam.to_str().unwrap());
    assert_eq!(
        ours, bt,
        "output diverged from bedtools\nours:\n{ours}\nbedtools:\n{bt}"
    );
}

#[test]
fn line_count_matches_bedtools() {
    if !bedtools_available() {
        eprintln!("skip: bedtools not found");
        return;
    }
    let bam = golden("small.bam");
    let ours = run_ours(bam.to_str().unwrap());
    let bt = run_bedtools(bam.to_str().unwrap());
    assert_eq!(
        ours.lines().count(),
        bt.lines().count(),
        "line count mismatch"
    );
}

#[test]
fn bed6_columns_correct() {
    let bam = golden("small.bam");
    let got = run_ours(bam.to_str().unwrap());
    let first = got.lines().next().expect("no output");
    let cols: Vec<&str> = first.split('\t').collect();
    assert_eq!(cols.len(), 6, "expected 6 BED columns, got {}", cols.len());
    // start < end
    let start: u64 = cols[1].parse().expect("start not numeric");
    let end: u64 = cols[2].parse().expect("end not numeric");
    assert!(start < end, "start ({start}) must be < end ({end})");
    // strand is + or -
    assert!(cols[5] == "+" || cols[5] == "-", "strand must be + or -");
}
