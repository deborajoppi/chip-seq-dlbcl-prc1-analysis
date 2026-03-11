#!/usr/bin/env Rscript

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1"
)

required <- c("fgsea", "msigdbr", "ggplot2")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  stop(
    sprintf("Missing required R packages for GSEA: %s", paste(missing, collapse = ", ")),
    call. = FALSE
  )
}

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

annotation_path <- "results/annotations/BCOR_KDM2B_overlap_no_H3K27me3.annotated.tsv"
if (!file.exists(annotation_path)) {
  stop("Missing annotation file for GSEA input. Run make annotate first.", call. = FALSE)
}

dat <- read.delim(annotation_path, stringsAsFactors = FALSE, check.names = FALSE)
genes <- trimws(dat$`Gene Name`)
genes[genes == "" | genes == "NA"] <- NA_character_

promoter_keep <- grepl("promoter|tss", tolower(dat$Annotation))
exclude_genes <- c("BCOR", "KDM2B")

all_counts <- sort(table(genes[!is.na(genes) & !(genes %in% exclude_genes)]), decreasing = TRUE)
promoter_counts <- sort(table(genes[promoter_keep & !is.na(genes) & !(genes %in% exclude_genes)]), decreasing = TRUE)

all_gene_names <- sort(unique(names(all_counts)))
rank_df <- data.frame(
  gene = all_gene_names,
  promoter_peak_count = 0,
  total_peak_count = 0,
  stringsAsFactors = FALSE
)
rank_df$promoter_peak_count[match(names(promoter_counts), rank_df$gene)] <- as.integer(promoter_counts)
rank_df$total_peak_count[match(names(all_counts), rank_df$gene)] <- as.integer(all_counts)
rank_df$rank_score <- rank_df$promoter_peak_count + (rank_df$total_peak_count / 1000)
rank_df <- rank_df[order(rank_df$rank_score, decreasing = TRUE, rank_df$gene), ]
rank_df$rank_score <- rank_df$rank_score + rev(seq_len(nrow(rank_df))) / 1e9

write.csv(rank_df, "results/tables/gsea_ranked_genes.csv", row.names = FALSE)

stats <- rank_df$rank_score
names(stats) <- rank_df$gene

hallmark <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H")
pathways <- split(hallmark$gene_symbol, hallmark$gs_name)

fgsea_res <- fgsea::fgsea(
  pathways = pathways,
  stats = stats,
  minSize = 10,
  maxSize = 500,
  scoreType = "pos"
)

fgsea_res <- fgsea_res[order(fgsea_res$padj, -abs(fgsea_res$NES)), ]
fgsea_out <- as.data.frame(fgsea_res)
if ("leadingEdge" %in% names(fgsea_out)) {
  fgsea_out$leadingEdge <- vapply(fgsea_out$leadingEdge, function(x) paste(x, collapse = ";"), character(1))
}
write.csv(fgsea_out, "results/tables/gsea_hallmark_promoter_overlap.csv", row.names = FALSE)

plot_df <- head(fgsea_out[is.finite(fgsea_out$NES) & !is.na(fgsea_out$padj), c("pathway", "NES", "padj")], 15)
plot_df$pathway <- sub("^HALLMARK_", "", plot_df$pathway)
plot_df$pathway <- factor(plot_df$pathway, levels = rev(plot_df$pathway))

p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = NES, y = pathway, fill = -log10(padj))) +
  ggplot2::geom_col(width = 0.75) +
  ggplot2::scale_fill_gradient(low = "#F4D35E", high = "#7F1D1D") +
  ggplot2::labs(
    title = "Hallmark GSEA on promoter-focused shared non-H3K27me3 targets",
    x = "Normalized enrichment score",
    y = NULL,
    fill = "-log10(adj. p)"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    panel.grid.major.y = ggplot2::element_blank(),
    legend.position = "right"
  )

ggplot2::ggsave(
  filename = "results/figures/gsea_hallmark_promoter_overlap.png",
  plot = p,
  width = 10,
  height = 7,
  dpi = 180
)
