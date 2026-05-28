use criterion::{Criterion, criterion_group, criterion_main};
use std::path::PathBuf;
use std::process::Command;

fn fixture() -> PathBuf {
    let env = std::env::var("BCMR_BENCH_DATA")
        .unwrap_or_else(|_| "/Volumes/Zane's HDD/rsomics-fixtures".to_string());
    PathBuf::from(env).join("medium.bam")
}

fn bench(c: &mut Criterion) {
    let bam = fixture();
    if !bam.exists() {
        eprintln!("skip: fixture not found at {}", bam.display());
        return;
    }

    let ours = env!("CARGO_BIN_EXE_rsomics-bam-to-bed");
    let mut group = c.benchmark_group("bam_to_bed");
    group.sample_size(10);

    group.bench_function("rsomics-bam-to-bed", |bm| {
        bm.iter(|| {
            let out = Command::new(ours)
                .args(["-i"])
                .arg(&bam)
                .output()
                .expect("ours run");
            assert!(
                out.status.success(),
                "{}",
                String::from_utf8_lossy(&out.stderr)
            );
        });
    });

    if Command::new("bedtools")
        .arg("--version")
        .output()
        .is_ok_and(|o| o.status.success())
    {
        group.bench_function("bedtools-bamtobed", |bm| {
            bm.iter(|| {
                let out = Command::new("bedtools")
                    .args(["bamtobed", "-i"])
                    .arg(&bam)
                    .output()
                    .expect("bedtools run");
                assert!(
                    out.status.success(),
                    "{}",
                    String::from_utf8_lossy(&out.stderr)
                );
            });
        });
    } else {
        eprintln!("bedtools not on PATH — skipping upstream comparison");
    }

    group.finish();
}

criterion_group!(benches, bench);
criterion_main!(benches);
