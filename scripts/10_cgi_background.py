#!/usr/bin/env python3
"""CpG-island background for the BCOR/KDM2B overlap.

KDM2B binds unmethylated CpG islands (CGIs) through its CXXC domain, so any
set defined by KDM2B binding is expected to be dominated by CGI promoters.
This script makes that expectation explicit, so downstream enrichment can be
tested against the right background:

1. Permutation test: what fraction of each peak set overlaps a CGI, compared
   with the same intervals shuffled across the genome (same chromosome,
   same lengths)?
2. Promoter universe: for every RefSeq gene, is its promoter (TSS +/- window)
   a CGI promoter, and is it bound by the shared BCOR-KDM2B regions (with and
   without the H3K27me3 filter)? scripts/05_run_enrichment.R uses this table
   to run over-representation analysis against an all-gene and a
   CGI-promoter-only background.

Outputs (results/tables/):
  cgi_overlap_permutation.csv
  cgi_promoter_summary.csv
  promoter_universe.tsv
"""

import argparse
import csv
import gzip
import random
import statistics
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ANN = ROOT / "data" / "raw" / "annotation"
PEAKS = ROOT / "results" / "tmp" / "normalized_peaks"
OVL = ROOT / "results" / "overlaps"
TABLES = ROOT / "results" / "tables"


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--chrominfo", default=ANN / "hg18_chromInfo.txt.gz", type=Path)
    p.add_argument("--cgi", default=ANN / "hg18_cpgIslandExt.txt.gz", type=Path)
    p.add_argument("--refgene", default=ANN / "hg18_refGene.txt.gz", type=Path)
    p.add_argument("--bcor", default=PEAKS / "BCOR_LY1_hg18.bed", type=Path)
    p.add_argument("--kdm2b", default=PEAKS / "KDM2B_LY1_hg18.bed", type=Path)
    p.add_argument("--shared", default=OVL / "bcor_kdm2b_overlap.bed", type=Path)
    p.add_argument("--shared-no-k27", default=OVL / "bcor_kdm2b_overlap_without_h3k27me3.bed", type=Path)
    p.add_argument("--outdir", default=TABLES, type=Path)
    p.add_argument("--n-shuffles", type=int, default=100)
    p.add_argument("--promoter-window", type=int, default=1000, help="bp either side of the TSS")
    p.add_argument("--seed", type=int, default=1)
    return p.parse_args()


def is_primary(chrom: str) -> bool:
    return chrom.startswith("chr") and "_" not in chrom and chrom != "chrM"


def open_text(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()


def ucsc_rows(path: Path):
    """Yield rows of a UCSC table dump, dropping the leading 'bin' column if present."""
    with open_text(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[0].isdigit() and len(f) > 1 and f[1].startswith("chr"):
                f = f[1:]
            elif f[0].isdigit() and len(f) > 2 and f[2].startswith("chr"):
                f = f[1:]  # refGene: bin, name, chrom, ...
            yield f


def write_genome(chrominfo: Path, out: Path) -> set:
    chroms = set()
    with open_text(chrominfo) as handle, out.open("w") as dst:
        for line in handle:
            chrom, size = line.split("\t")[:2]
            if is_primary(chrom):
                dst.write(f"{chrom}\t{int(size)}\n")
                chroms.add(chrom)
    return chroms


def write_cgi(cgi: Path, out: Path, chroms: set) -> int:
    rows = []
    for f in ucsc_rows(cgi):
        chrom, start, end = f[0], int(f[1]), int(f[2])
        if chrom in chroms:
            rows.append((chrom, start, end))
    rows.sort()
    with out.open("w") as dst:
        for r in rows:
            dst.write("%s\t%d\t%d\n" % r)
    return len(rows)


def write_promoters(refgene: Path, out: Path, chroms: set, window: int) -> None:
    """One line per transcript TSS window: chrom, start, end, refseq_id, symbol."""
    rows = []
    for f in ucsc_rows(refgene):
        # name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount,
        # exonStarts, exonEnds, score, name2, ...
        refseq, chrom, strand = f[0], f[1], f[2]
        tx_start, tx_end, symbol = int(f[3]), int(f[4]), f[11]
        if chrom not in chroms:
            continue
        tss = tx_start if strand == "+" else tx_end
        rows.append((chrom, max(0, tss - window), tss + window, refseq.split(".")[0], symbol))
    rows.sort()
    with out.open("w") as dst:
        for r in rows:
            dst.write("%s\t%d\t%d\t%s\t%s\n" % r)


def bed3_primary(src: Path, out: Path, chroms: set) -> int:
    n = 0
    with src.open() as handle, out.open("w") as dst:
        for line in handle:
            f = line.split("\t")
            if len(f) >= 3 and f[0] in chroms:
                dst.write(f"{f[0]}\t{f[1]}\t{f[2].strip()}\n")
                n += 1
    return n


def bedtools(*args) -> str:
    return subprocess.run(["bedtools", *map(str, args)], check=True, capture_output=True, text=True).stdout


def frac_overlapping(bed: Path, cgi: Path) -> float:
    total = sum(1 for _ in bed.open())
    if total == 0:
        return 0.0
    hit = sum(1 for line in bedtools("intersect", "-u", "-a", bed, "-b", cgi).splitlines() if line.strip())
    return hit / total


def permutation(name: str, bed: Path, cgi: Path, genome: Path, n: int, rng: random.Random, tmp: Path) -> dict:
    observed = frac_overlapping(bed, cgi)
    shuffled_path = tmp / "shuffled.bed"
    null = []
    for _ in range(n):
        shuffled_path.write_text(
            bedtools("shuffle", "-i", bed, "-g", genome, "-chrom", "-seed", rng.randint(1, 2**31 - 1))
        )
        null.append(frac_overlapping(shuffled_path, cgi))
    mean_null = statistics.fmean(null)
    return {
        "peak_set": name,
        "n_intervals": sum(1 for _ in bed.open()),
        "observed_fraction_on_cgi": round(observed, 4),
        "shuffled_mean_fraction": round(mean_null, 4),
        "shuffled_sd": round(statistics.pstdev(null), 4),
        "fold_enrichment": round(observed / mean_null, 2) if mean_null > 0 else "",
        "empirical_p": round((sum(x >= observed for x in null) + 1) / (n + 1), 4),
        "n_shuffles": n,
    }


def genes_hit(promoters: Path, regions: Path) -> set:
    out = bedtools("intersect", "-u", "-a", promoters, "-b", regions)
    return {line.split("\t")[4].strip() for line in out.splitlines() if line.strip()}


def main() -> None:
    args = parse_args()
    for path in (args.chrominfo, args.cgi, args.refgene, args.bcor, args.kdm2b, args.shared, args.shared_no_k27):
        if not Path(path).exists():
            raise SystemExit(f"Missing input: {path}. Run make fetch / prepare / overlap first.")
    args.outdir.mkdir(parents=True, exist_ok=True)
    rng = random.Random(args.seed)

    with tempfile.TemporaryDirectory(prefix="cgi-bg-") as tmpdir:
        tmp = Path(tmpdir)
        genome = tmp / "genome.txt"
        chroms = write_genome(args.chrominfo, genome)
        cgi = tmp / "cgi.bed"
        n_cgi = write_cgi(args.cgi, cgi, chroms)
        promoters = tmp / "promoters.bed"
        write_promoters(args.refgene, promoters, chroms, args.promoter_window)
        print(f"Loaded {n_cgi} CpG islands on {len(chroms)} primary chromosomes")

        sets = {
            "BCOR peaks": args.bcor,
            "KDM2B peaks": args.kdm2b,
            "Shared BCOR-KDM2B regions": args.shared,
            "Shared regions outside H3K27me3 mask": args.shared_no_k27,
        }
        clean = {}
        for i, (name, src) in enumerate(sets.items()):
            clean[name] = tmp / f"set{i}.bed"
            bed3_primary(src, clean[name], chroms)

        # 1. Permutation test
        perm_rows = [permutation(name, bed, cgi, genome, args.n_shuffles, rng, tmp) for name, bed in clean.items()]
        with (args.outdir / "cgi_overlap_permutation.csv").open("w", newline="") as handle:
            w = csv.DictWriter(handle, fieldnames=list(perm_rows[0].keys()), lineterminator="\n")
            w.writeheader()
            w.writerows(perm_rows)

        # 2. Gene-level promoter universe
        refseqs = defaultdict(set)
        with promoters.open() as handle:
            for line in handle:
                f = line.rstrip("\n").split("\t")
                refseqs[f[4]].add(f[3])
        cgi_genes = genes_hit(promoters, cgi)
        shared_genes = genes_hit(promoters, clean["Shared BCOR-KDM2B regions"])
        no_k27_genes = genes_hit(promoters, clean["Shared regions outside H3K27me3 mask"])

        with (args.outdir / "promoter_universe.tsv").open("w", newline="") as handle:
            w = csv.writer(handle, delimiter="\t", lineterminator="\n")
            w.writerow(["symbol_hg18", "refseq_ids", "cgi_promoter", "bound_shared", "bound_shared_no_h3k27me3"])
            for sym in sorted(refseqs):
                w.writerow([sym, ";".join(sorted(refseqs[sym])),
                            int(sym in cgi_genes), int(sym in shared_genes), int(sym in no_k27_genes)])

        all_genes = set(refseqs)
        summary = []
        for label, genes in [("All RefSeq genes", all_genes),
                             ("Promoter-bound by shared regions", shared_genes),
                             ("Promoter-bound by shared regions outside H3K27me3", no_k27_genes)]:
            n = len(genes)
            k = len(genes & cgi_genes)
            summary.append({"gene_set": label, "n_genes": n, "n_cgi_promoter": k,
                            "fraction_cgi_promoter": round(k / n, 4) if n else ""})
        with (args.outdir / "cgi_promoter_summary.csv").open("w", newline="") as handle:
            w = csv.DictWriter(handle, fieldnames=list(summary[0].keys()), lineterminator="\n")
            w.writeheader()
            w.writerows(summary)

    for r in perm_rows:
        print(f"{r['peak_set']}: {r['observed_fraction_on_cgi']:.1%} on CGI vs "
              f"{r['shuffled_mean_fraction']:.1%} shuffled ({r['fold_enrichment']}x, p={r['empirical_p']})")
    for r in summary:
        print(f"{r['gene_set']}: {r['n_cgi_promoter']}/{r['n_genes']} CGI promoters")
    print(f"Wrote tables to {args.outdir}/")


if __name__ == "__main__":
    main()
