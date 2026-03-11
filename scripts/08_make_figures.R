#!/usr/bin/env Rscript

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

classify_region <- function(x) {
  x <- tolower(x)
  if (grepl("promoter|tss", x)) return("promoter-TSS")
  if (grepl("intron", x)) return("intron")
  if (grepl("exon", x)) return("exon")
  if (grepl("intergenic", x)) return("intergenic")
  if (grepl("tts", x)) return("TTS")
  "other"
}

read_region_counts <- function(path) {
  dat <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  counts <- table(vapply(dat$Annotation, classify_region, character(1)))
  categories <- c("promoter-TSS", "intron", "exon", "intergenic", "TTS")
  out <- setNames(integer(length(categories)), categories)
  out[names(counts)] <- as.integer(counts)
  out
}

read_gene_counts <- function(path) {
  dat <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  genes <- trimws(dat$`Gene Name`)
  genes <- genes[genes != "" & genes != "NA"]
  sort(table(genes), decreasing = TRUE)
}

read_promoter_gene_counts <- function(path, exclude_genes = character()) {
  dat <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  keep <- grepl("promoter|tss", tolower(dat$Annotation))
  genes <- trimws(dat$`Gene Name`[keep])
  genes <- genes[genes != "" & genes != "NA" & !(genes %in% exclude_genes)]
  sort(table(genes), decreasing = TRUE)
}

read_set_gene_counts <- function(path, promoter_only = FALSE, exclude_genes = character()) {
  dat <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (promoter_only) {
    dat <- dat[grepl("promoter|tss", tolower(dat$Annotation)), , drop = FALSE]
  }
  genes <- trimws(dat$`Gene Name`)
  genes <- genes[genes != "" & genes != "NA" & !(genes %in% exclude_genes)]
  sort(table(genes), decreasing = TRUE)
}

make_venn <- function() {
  counts <- read.csv("results/tables/thesis_count_comparison.csv", stringsAsFactors = FALSE)
  bcor_unique <- counts$rerun_observed[counts$metric == "BCOR unique"]
  kdm2b_unique <- counts$rerun_observed[counts$metric == "KDM2B unique"]
  overlap <- counts$rerun_observed[counts$metric == "BCOR and KDM2B overlap"]

  png("results/figures/figure31a_venn_rerun.png", width = 1600, height = 1200, res = 180)
  par(mar = c(2, 2, 4, 2))
  plot.new()
  plot.window(xlim = c(0, 10), ylim = c(0, 6))
  symbols(4, 3, circles = 2.2, inches = FALSE, add = TRUE, bg = rgb(0.85, 0.35, 0.25, 0.35), fg = "#B33A2B")
  symbols(6, 3, circles = 2.2, inches = FALSE, add = TRUE, bg = rgb(0.20, 0.45, 0.80, 0.35), fg = "#2F5AA8")
  text(2.8, 5.4, "BCOR", cex = 1.8, font = 2, col = "#7C2118")
  text(7.2, 5.4, "KDM2B", cex = 1.8, font = 2, col = "#17376C")
  text(2.6, 3.0, format(bcor_unique, big.mark = ","), cex = 2.2, font = 2)
  text(5.0, 3.0, format(overlap, big.mark = ","), cex = 2.2, font = 2)
  text(7.4, 3.0, format(kdm2b_unique, big.mark = ","), cex = 2.2, font = 2)
  title("BCOR and KDM2B overlap in the public rerun", cex.main = 1.8, font.main = 2)
  mtext("Shared regions are distinct BCOR-KDM2B genomic intersections", side = 3, line = 0.5, cex = 1.0)
  dev.off()
}

make_region_plot <- function() {
  counts <- rbind(
    BCOR = read_region_counts("results/annotations/BCOR.annotated.tsv"),
    KDM2B = read_region_counts("results/annotations/KDM2B.annotated.tsv"),
    Overlap = read_region_counts("results/annotations/BCOR_KDM2B_overlap.annotated.tsv")
  )
  counts <- t(counts[, c("promoter-TSS", "intron", "exon", "intergenic", "TTS")])
  cols <- c("#C8523B", "#D98B3A", "#E8C547", "#5B8E7D", "#3E6FB6")

  png("results/figures/figure31b_genomic_regions_rerun.png", width = 1800, height = 1200, res = 180)
  par(mar = c(6, 5, 4, 2))
  barplot(
    counts,
    beside = FALSE,
    col = cols,
    border = NA,
    las = 1,
    ylab = "Peak count",
    main = "Genomic distribution of BCOR, KDM2B, and shared regions"
  )
  legend("topright", legend = rownames(counts), fill = cols, bty = "n", cex = 1.0)
  dev.off()
}

make_threshold_plot <- function() {
  if (!file.exists("results/tables/h3k27me3_threshold_sweep.csv")) return(invisible(NULL))
  dat <- read.csv("results/tables/h3k27me3_threshold_sweep.csv", stringsAsFactors = FALSE)

  png("results/figures/h3k27me3_threshold_sensitivity.png", width = 1800, height = 1200, res = 180)
  par(mar = c(5, 5, 4, 2))
  plot(
    dat$min_signal, dat$overlap_without_h3k27me3,
    type = "b", pch = 19, lwd = 3, col = "#2F5AA8",
    xlab = "Minimum H3K27me3 signal threshold",
    ylab = "Overlap regions retained",
    ylim = range(c(dat$overlap_without_h3k27me3, dat$overlap_with_h3k27me3)),
    main = "Sensitivity of PRC2 exclusion to H3K27me3 threshold"
  )
  lines(dat$min_signal, dat$overlap_with_h3k27me3, type = "b", pch = 17, lwd = 3, col = "#B33A2B")
  legend(
    "topright",
    legend = c("Retained after exclusion", "Excluded by H3K27me3"),
    col = c("#2F5AA8", "#B33A2B"),
    lwd = 3,
    pch = c(19, 17),
    bty = "n"
  )
  dev.off()
}

make_top_target_barplot <- function() {
  counts <- read_gene_counts("results/annotations/BCOR_KDM2B_overlap_no_H3K27me3.annotated.tsv")
  top_counts <- head(counts, 20)
  top_df <- data.frame(
    gene = names(top_counts),
    peak_count = as.integer(top_counts),
    stringsAsFactors = FALSE
  )
  write.csv(top_df, "results/tables/top_candidate_gene_peak_counts.csv", row.names = FALSE)

  png("results/figures/top_candidate_genes_barplot.png", width = 1800, height = 1400, res = 180)
  par(mar = c(6, 10, 4, 2))
  barplot(
    rev(top_df$peak_count),
    horiz = TRUE,
    names.arg = rev(top_df$gene),
    las = 1,
    col = "#C8523B",
    border = NA,
    xlab = "Number of assigned shared non-H3K27me3 peaks",
    main = "Top candidate target genes by retained shared peak count"
  )
  dev.off()
}

make_gene_set_heatmap <- function() {
  set_paths <- c(
    BCOR = "results/annotations/BCOR.annotated.tsv",
    KDM2B = "results/annotations/KDM2B.annotated.tsv",
    Overlap = "results/annotations/BCOR_KDM2B_overlap.annotated.tsv",
    Overlap_no_H3K27me3 = "results/annotations/BCOR_KDM2B_overlap_no_H3K27me3.annotated.tsv"
  )

  top_counts <- read_promoter_gene_counts(
    "results/annotations/BCOR_KDM2B_overlap_no_H3K27me3.annotated.tsv",
    exclude_genes = c("BCOR", "KDM2B")
  )
  top_genes <- names(head(top_counts, 20))

  mat <- matrix(0, nrow = length(top_genes), ncol = length(set_paths), dimnames = list(top_genes, names(set_paths)))
  for (set_name in names(set_paths)) {
    counts <- read_set_gene_counts(set_paths[[set_name]], promoter_only = TRUE, exclude_genes = c("BCOR", "KDM2B"))
    mat[, set_name] <- as.integer(counts[top_genes])
    mat[is.na(mat[, set_name]), set_name] <- 0L
  }

  write.csv(
    data.frame(gene = rownames(mat), mat, row.names = NULL, check.names = FALSE),
    "results/tables/candidate_gene_set_matrix.csv",
    row.names = FALSE
  )

  z <- log1p(mat)
  cols <- colorRampPalette(c("#FFF8E7", "#F4D35E", "#EE964B", "#D1495B", "#7F1D1D"))(100)

  png("results/figures/top_candidate_genes_heatmap.png", width = 1600, height = 1800, res = 180)
  par(mar = c(10, 12, 4, 2))
  image(
    x = seq_len(ncol(z)),
    y = seq_len(nrow(z)),
    z = t(z[nrow(z):1, , drop = FALSE]),
    col = cols,
    axes = FALSE,
    xlab = "",
    ylab = "",
    main = "Top promoter-bound candidate genes across peak sets"
  )
  axis(1, at = seq_len(ncol(z)), labels = colnames(z), las = 2, cex.axis = 0.9)
  axis(2, at = seq_len(nrow(z)), labels = rev(rownames(z)), las = 2, cex.axis = 0.75)
  for (i in seq_len(nrow(mat))) {
    for (j in seq_len(ncol(mat))) {
      text(j, nrow(mat) - i + 1, labels = mat[i, j], cex = 0.65, col = ifelse(z[i, j] > median(z), "white", "#3A2A16"))
    }
  }
  box()
  dev.off()
}

make_promoter_target_barplot <- function() {
  counts <- read_promoter_gene_counts(
    "results/annotations/BCOR_KDM2B_overlap_no_H3K27me3.annotated.tsv",
    exclude_genes = c("BCOR", "KDM2B")
  )
  top_counts <- head(counts, 20)
  top_df <- data.frame(
    gene = names(top_counts),
    promoter_peak_count = as.integer(top_counts),
    stringsAsFactors = FALSE
  )
  write.csv(top_df, "results/tables/top_promoter_candidate_gene_peak_counts.csv", row.names = FALSE)

  png("results/figures/top_promoter_candidate_genes_barplot.png", width = 1800, height = 1400, res = 180)
  par(mar = c(6, 10, 4, 2))
  barplot(
    rev(top_df$promoter_peak_count),
    horiz = TRUE,
    names.arg = rev(top_df$gene),
    las = 1,
    col = "#2F5AA8",
    border = NA,
    xlab = "Number of promoter-assigned shared non-H3K27me3 peaks",
    main = "Top promoter-bound candidate genes excluding BCOR and KDM2B"
  )
  dev.off()
}

make_venn()
make_region_plot()
make_threshold_plot()
make_top_target_barplot()
make_gene_set_heatmap()
make_promoter_target_barplot()
