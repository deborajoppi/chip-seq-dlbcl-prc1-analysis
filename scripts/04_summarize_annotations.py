#!/usr/bin/env python3

import csv
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OVERLAP_DIR = ROOT / "results" / "overlaps"
ANNOTATION_DIR = ROOT / "results" / "annotations"
TABLE_DIR = ROOT / "results" / "tables"


def count_lines(path: Path) -> int:
    if not path.exists():
        return 0
    with path.open() as handle:
        return sum(1 for line in handle if line.strip())


def classify_annotation(value: str) -> str:
    text = value.lower()
    if "promoter" in text or "tss" in text:
        return "promoter-TSS"
    if "intron" in text:
        return "intron"
    if "exon" in text:
        return "exon"
    if "tts" in text:
        return "TTS"
    if "intergenic" in text:
        return "intergenic"
    return "other"


def write_peak_totals() -> None:
    rows = [
        ("BCOR", count_lines(ROOT / "results" / "tmp" / "normalized_peaks" / "BCOR_LY1_hg18.bed")),
        ("KDM2B", count_lines(ROOT / "results" / "tmp" / "normalized_peaks" / "KDM2B_LY1_hg18.bed")),
        ("BCOR and KDM2B overlap", count_lines(OVERLAP_DIR / "bcor_kdm2b_overlap.bed")),
        ("BCOR unique vs KDM2B", count_lines(OVERLAP_DIR / "bcor_unique_vs_kdm2b.bed")),
        ("KDM2B unique vs BCOR", count_lines(OVERLAP_DIR / "kdm2b_unique_vs_bcor.bed")),
        ("BCOR and KDM2B overlap with H3K27me3", count_lines(OVERLAP_DIR / "bcor_kdm2b_overlap_with_h3k27me3.bed")),
        ("BCOR and KDM2B overlap without H3K27me3", count_lines(OVERLAP_DIR / "bcor_kdm2b_overlap_without_h3k27me3.bed")),
    ]

    with (TABLE_DIR / "peak_totals.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["target_set", "total_peaks"])
        writer.writerows(rows)


def read_annotation_table(path: Path):
    if not path.exists():
        raise FileNotFoundError(f"Missing annotation file: {path}")

    with path.open() as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, None)
        if header is None:
            return [], None, None

        annotation_idx = None
        gene_idx = None
        for i, name in enumerate(header):
            lower = name.lower()
            if lower == "annotation":
                annotation_idx = i
            if lower in {"gene name", "nearest promoterid", "nearest refseq", "ensembl", "alias"} and gene_idx is None:
                gene_idx = i

        rows = list(reader)
        return rows, annotation_idx, gene_idx


def write_region_summary() -> None:
    rows, annotation_idx, _ = read_annotation_table(ANNOTATION_DIR / "BCOR_KDM2B_overlap.annotated.tsv")
    counts = Counter()
    if annotation_idx is not None:
        for row in rows:
            if annotation_idx < len(row):
                counts[classify_annotation(row[annotation_idx])] += 1

    with (TABLE_DIR / "genomic_region_summary.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["peak_set", "region_class", "count"])
        for key in ["promoter-TSS", "intron", "exon", "intergenic", "TTS", "other"]:
            writer.writerow(["BCOR and KDM2B overlap", key, counts.get(key, 0)])


def write_candidate_genes() -> None:
    rows, _, gene_idx = read_annotation_table(ANNOTATION_DIR / "BCOR_KDM2B_overlap_no_H3K27me3.annotated.tsv")
    genes = []
    if gene_idx is not None:
        for row in rows:
            if gene_idx < len(row):
                gene = row[gene_idx].strip()
                if gene and gene != "NA":
                    genes.append(gene)

    unique_genes = sorted(set(genes))
    with (TABLE_DIR / "candidate_target_genes.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["gene", "source"])
        for gene in unique_genes:
            writer.writerow([gene, "BCOR_KDM2B_overlap_without_H3K27me3"])


def main() -> None:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    write_peak_totals()
    write_region_summary()
    write_candidate_genes()
    print("Wrote summary tables to results/tables/")


if __name__ == "__main__":
    main()
