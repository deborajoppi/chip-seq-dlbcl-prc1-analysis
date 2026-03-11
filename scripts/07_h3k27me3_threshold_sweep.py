#!/usr/bin/env python3

import csv
import subprocess
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
WIG = ROOT / "data" / "raw" / "peaks" / "H3K27me3_LY1_hg18.source.wig.gz"
OVERLAP = ROOT / "results" / "overlaps" / "bcor_kdm2b_overlap.bed"
CURRENT_MASK = ROOT / "data" / "raw" / "peaks" / "H3K27me3_LY1_hg18.bed"
OUT = ROOT / "results" / "tables" / "h3k27me3_threshold_sweep.csv"
THRESHOLDS = [0.1, 0.5, 1.0]


def line_count(path: Path) -> int:
    with path.open() as handle:
        return sum(1 for _ in handle)


def intersect_count(mask: Path, keep_flag: str) -> int:
    cmd = ["bedtools", "intersect", "-a", str(OVERLAP), "-b", str(mask), keep_flag]
    proc = subprocess.run(cmd, check=True, capture_output=True, text=True)
    return sum(1 for line in proc.stdout.splitlines() if line.strip())


def build_mask(threshold: float, output: Path) -> Path:
    if threshold == 0.1 and CURRENT_MASK.exists():
        return CURRENT_MASK
    cached = Path("/tmp") / f"h3k27me3_{str(threshold).replace('.', 'p')}.bed"
    if cached.exists():
        return cached

    subprocess.run(
        [
            "python3",
            "scripts/wig_to_bed.py",
            "--input",
            str(WIG),
            "--output",
            str(output),
            "--min-signal",
            str(threshold),
        ],
        check=True,
        cwd=ROOT,
    )
    return output


def main() -> None:
    if not WIG.exists():
        raise SystemExit(f"Missing WIG source file: {WIG}")
    if not OVERLAP.exists():
        raise SystemExit(f"Missing overlap BED: {OVERLAP}. Run make overlap first.")

    rows = []
    with tempfile.TemporaryDirectory(prefix="h3k27me3-threshold-", dir="/tmp") as tmpdir:
        tmpdir = Path(tmpdir)
        for threshold in THRESHOLDS:
            mask_path = build_mask(threshold, tmpdir / f"h3k27me3_{str(threshold).replace('.', 'p')}.bed")
            rows.append(
                {
                    "min_signal": threshold,
                    "mask_intervals": line_count(mask_path),
                    "overlap_with_h3k27me3": intersect_count(mask_path, "-u"),
                    "overlap_without_h3k27me3": intersect_count(mask_path, "-v"),
                }
            )

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["min_signal", "mask_intervals", "overlap_with_h3k27me3", "overlap_without_h3k27me3"],
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
