#!/usr/bin/env Rscript

required <- c("clusterProfiler", "org.Hs.eg.db")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  stop(
    sprintf(
      "Missing required R packages for GO enrichment: %s",
      paste(missing, collapse = ", ")
    ),
    call. = FALSE
  )
}

genes <- read.csv("results/tables/candidate_target_genes.csv", stringsAsFactors = FALSE)
if (nrow(genes) == 0) {
  stop("No candidate genes found in results/tables/candidate_target_genes.csv", call. = FALSE)
}

ego <- clusterProfiler::enrichGO(
  gene = unique(genes$gene),
  OrgDb = org.Hs.eg.db::org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  readable = TRUE
)

out <- as.data.frame(ego)
write.csv(out, "results/tables/go_enrichment.csv", row.names = FALSE)
message("Wrote results/tables/go_enrichment.csv")
