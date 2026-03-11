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
    bcor_total = count_lines(ROOT / "results" / "tmp" / "normalized_peaks" / "BCOR_LY1_hg18.bed")
    kdm2b_total = count_lines(ROOT / "results" / "tmp" / "normalized_peaks" / "KDM2B_LY1_hg18.bed")
    overlap_total = count_lines(OVERLAP_DIR / "bcor_kdm2b_overlap.bed")

    rows = [
        ("BCOR", bcor_total),
        ("KDM2B", kdm2b_total),
        ("BCOR and KDM2B overlap", overlap_total),
        ("BCOR unique (thesis-style: total minus overlap)", bcor_total - overlap_total),
        ("KDM2B unique (thesis-style: total minus overlap)", kdm2b_total - overlap_total),
        ("BCOR non-overlap strict", count_lines(OVERLAP_DIR / "bcor_nonoverlap_strict.bed")),
        ("KDM2B non-overlap strict", count_lines(OVERLAP_DIR / "kdm2b_nonoverlap_strict.bed")),
        ("BCOR and KDM2B overlap with H3K27me3", count_lines(OVERLAP_DIR / "bcor_kdm2b_overlap_with_h3k27me3.bed")),
        ("BCOR and KDM2B overlap without H3K27me3", count_lines(OVERLAP_DIR / "bcor_kdm2b_overlap_without_h3k27me3.bed")),
    ]

    with (TABLE_DIR / "peak_totals.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["target_set", "total_peaks"])
        writer.writerows(rows)

    with (TABLE_DIR / "thesis_count_comparison.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["metric", "thesis_reported", "rerun_observed"])
        writer.writerow(["BCOR", 17549, bcor_total])
        writer.writerow(["KDM2B", 26471, kdm2b_total])
        writer.writerow(["BCOR and KDM2B overlap", 14546, overlap_total])
        writer.writerow(["BCOR unique", 3002, bcor_total - overlap_total])
        writer.writerow(["KDM2B unique", 11925, kdm2b_total - overlap_total])


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
        preferred_gene_columns = [
            "gene name",
            "gene alias",
            "nearest promoterid",
            "nearest refseq",
            "nearest ensembl",
            "entrez id",
            "alias",
        ]
        header_index = {name.lower(): i for i, name in enumerate(header)}

        for i, name in enumerate(header):
            lower = name.lower()
            if lower == "annotation":
                annotation_idx = i

        for column_name in preferred_gene_columns:
            if column_name in header_index:
                gene_idx = header_index[column_name]
                break

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
