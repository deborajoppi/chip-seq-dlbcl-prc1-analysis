#!/usr/bin/env Rscript
# Over-representation analysis (GO BP + MSigDB Hallmark) for genes whose
# promoter is bound by shared BCOR-KDM2B regions outside the H3K27me3 mask.
#
# Each test is run against two explicit backgrounds:
#   all_genes     - every RefSeq gene with a promoter on a primary chromosome
#   cgi_promoters - only genes with a CpG-island promoter (targets restricted
#                   to CGI promoters as well)
# Terms that are enriched against all_genes but not against cgi_promoters are
# most likely explained by KDM2B's general preference for CpG islands rather
# than by PRC1.1-specific targeting.
#
# Input:  results/tables/promoter_universe.tsv (from scripts/10_cgi_background.py)
# Output: results/tables/enrichment_<db>_<background>.csv
#         results/tables/enrichment_summary.csv

required <- c("clusterProfiler", "org.Hs.eg.db", "msigdbr", "AnnotationDbi")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  stop(sprintf("Missing R packages: %s (see environment.yml)", paste(missing, collapse = ", ")), call. = FALSE)
}

universe_path <- "results/tables/promoter_universe.tsv"
if (!file.exists(universe_path)) {
  stop("Missing promoter_universe.tsv. Run make cgi first.", call. = FALSE)
}

set.seed(1)
u <- read.delim(universe_path, stringsAsFactors = FALSE)

# hg18-era gene symbols are often outdated, so map through RefSeq accessions.
refseq_long <- data.frame(
  symbol_hg18 = rep(u$symbol_hg18, lengths(strsplit(u$refseq_ids, ";"))),
  refseq = unlist(strsplit(u$refseq_ids, ";")),
  stringsAsFactors = FALSE
)
map <- AnnotationDbi::select(
  org.Hs.eg.db::org.Hs.eg.db, keys = unique(refseq_long$refseq),
  keytype = "REFSEQ", columns = "ENTREZID"
)
map <- map[!is.na(map$ENTREZID), ]
refseq_long$entrez <- map$ENTREZID[match(refseq_long$refseq, map$REFSEQ)]
gene_entrez <- tapply(refseq_long$entrez, refseq_long$symbol_hg18, function(x) {
  x <- x[!is.na(x)]
  if (length(x)) names(sort(table(x), decreasing = TRUE))[1] else NA_character_
})
u$entrez <- gene_entrez[u$symbol_hg18]
message(sprintf("Mapped %d of %d hg18 genes to Entrez IDs", sum(!is.na(u$entrez)), nrow(u)))
u <- u[!is.na(u$entrez) & !(u$symbol_hg18 %in% c("BCOR", "KDM2B", "FBXL10")), ]

backgrounds <- list(
  all_genes = list(
    universe = unique(u$entrez),
    targets = unique(u$entrez[u$bound_shared_no_h3k27me3 == 1])
  ),
  cgi_promoters = list(
    universe = unique(u$entrez[u$cgi_promoter == 1]),
    targets = unique(u$entrez[u$bound_shared_no_h3k27me3 == 1 & u$cgi_promoter == 1])
  )
)

hallmark <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H")
hallmark_t2g <- unique(hallmark[, c("gs_name", "ncbi_gene")])
names(hallmark_t2g) <- c("term", "gene")
hallmark_t2g$gene <- as.character(hallmark_t2g$gene)

summary_rows <- list()
for (bg in names(backgrounds)) {
  b <- backgrounds[[bg]]
  message(sprintf("[%s] %d targets in a universe of %d genes", bg, length(b$targets), length(b$universe)))

  go <- clusterProfiler::enrichGO(
    gene = b$targets, universe = b$universe, OrgDb = org.Hs.eg.db::org.Hs.eg.db,
    keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH",
    pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE
  )
  hm <- clusterProfiler::enricher(
    gene = b$targets, universe = b$universe, TERM2GENE = hallmark_t2g,
    pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1
  )

  for (res in list(list(db = "go_bp", obj = go), list(db = "hallmark", obj = hm))) {
    df <- if (is.null(res$obj)) data.frame() else as.data.frame(res$obj)
    write.csv(df, sprintf("results/tables/enrichment_%s_%s.csv", res$db, bg), row.names = FALSE)
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      database = res$db, background = bg,
      n_targets = length(b$targets), n_universe = length(b$universe),
      n_terms_padj_lt_0.05 = if (nrow(df)) sum(df$p.adjust < 0.05) else 0,
      top_terms = if (nrow(df)) paste(head(df$Description[df$p.adjust < 0.05], 5), collapse = "; ") else ""
    )
  }
}

summary_df <- do.call(rbind, summary_rows)
write.csv(summary_df, "results/tables/enrichment_summary.csv", row.names = FALSE)
print(summary_df[, 1:5])
message("Wrote results/tables/enrichment_*.csv")
