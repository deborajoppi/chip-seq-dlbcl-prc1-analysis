#!/usr/bin/env python3
"""Sensitivity of the H3K27me3 exclusion step to the WIG signal threshold.

For each threshold the mask is rebuilt from the source WIG (no caching
between runs), and we report how many shared BCOR-KDM2B regions are kept
or excluded, plus how much of the genome the mask covers. A mask covering a
large share of the genome means "outside H3K27me3" mostly reflects absence
of reads rather than absence of the mark.
"""

import argparse
import csv
import gzip
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
WIG = ROOT / "data" / "raw" / "peaks" / "H3K27me3_LY1_hg18.source.wig.gz"
CHROMINFO = ROOT / "data" / "raw" / "annotation" / "hg18_chromInfo.txt.gz"
OVERLAP = ROOT / "results" / "overlaps" / "bcor_kdm2b_overlap.bed"
OUT = ROOT / "results" / "tables" / "h3k27me3_threshold_sweep.csv"
DEFAULT_THRESHOLDS = [0.1, 0.25, 0.5, 1.0, 2.0]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--thresholds", type=float, nargs="+", default=DEFAULT_THRESHOLDS)
    return p.parse_args()


def genome_size(path: Path):
    """Sum of primary chromosome lengths (excludes _random, chrM, haplotypes)."""
    if not path.exists():
        return None
    total = 0
    with gzip.open(path, "rt") as handle:
        for line in handle:
            chrom, size = line.split("\t")[:2]
            if "_" in chrom or chrom == "chrM":
                continue
            total += int(size)
    return total


def mask_stats(mask: Path):
    intervals = 0
    bp = 0
    with mask.open() as handle:
        for line in handle:
            if not line.strip():
                continue
            f = line.split("\t")
            intervals += 1
            bp += int(f[2]) - int(f[1])
    return intervals, bp


def intersect_count(mask: Path, flag: str) -> int:
    cmd = ["bedtools", "intersect", "-a", str(OVERLAP), "-b", str(mask), flag]
    proc = subprocess.run(cmd, check=True, capture_output=True, text=True)
    return sum(1 for line in proc.stdout.splitlines() if line.strip())


def main() -> None:
    args = parse_args()
    if not WIG.exists():
        raise SystemExit(f"Missing WIG source file: {WIG}. Run make fetch first.")
    if not OVERLAP.exists():
        raise SystemExit(f"Missing overlap BED: {OVERLAP}. Run make overlap first.")

    gsize = genome_size(CHROMINFO)
    if gsize is None:
        print(f"Note: {CHROMINFO} not found; genome fraction will be left blank.")

    rows = []
    with tempfile.TemporaryDirectory(prefix="h3k27me3-sweep-") as tmp:
        for t in args.thresholds:
            mask = Path(tmp) / f"mask_{t}.bed"
            subprocess.run(
                ["python3", str(ROOT / "scripts" / "wig_to_bed.py"),
                 "--input", str(WIG), "--output", str(mask), "--min-signal", str(t)],
                check=True,
            )
            n_int, bp = mask_stats(mask)
            rows.append({
                "min_signal": t,
                "mask_intervals": n_int,
                "mask_bp": bp,
                "genome_fraction_masked": round(bp / gsize, 4) if gsize else "",
                "overlap_with_h3k27me3": intersect_count(mask, "-u"),
                "overlap_without_h3k27me3": intersect_count(mask, "-v"),
            })

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
