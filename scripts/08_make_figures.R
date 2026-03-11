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

make_venn()
make_region_plot()
make_threshold_plot()
