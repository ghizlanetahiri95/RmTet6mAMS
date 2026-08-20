# GO Biological Process and KEGG pathway enrichment analysis
#
# Identifies GO and KEGG terms enriched among genes commonly upregulated
# in at least 3 of 4 organs in untreated infected mice vs non-infected controls.
#
# Requires: vst_matrix and metadata in the environment (from 02_counts_to_vst.R)

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(AnnotationDbi)
  library(ggplot2)
  library(dplyr)
})

# -----------------------------------------------------------------------------
# Identify commonly upregulated genes (>= 3 organs)
# -----------------------------------------------------------------------------

organs <- c("Liver", "Lung", "Spleen", "Brain")
fc_threshold  <- 0.5   # minimum mean FC (VST scale) to call a gene upregulated
min_organs    <- 3     # gene must be upregulated in at least this many organs

up_by_organ <- lapply(organs, function(org) {
  s_notra <- metadata$Treatment == "Notra" & metadata$Organ == org
  s_noinf <- metadata$Treatment == "Noinf" & metadata$Organ == org
  fc <- rowMeans(vst_symbol[, s_notra, drop = FALSE]) -
        rowMeans(vst_symbol[, s_noinf, drop = FALSE])
  names(fc[fc > fc_threshold])
})
names(up_by_organ) <- organs

# Genes upregulated in >= min_organs
gene_counts   <- table(unlist(up_by_organ))
common_up     <- names(gene_counts[gene_counts >= min_organs])
cat("Genes upregulated in >=", min_organs, "organs:", length(common_up), "\n\n")

# Convert gene symbols to Entrez IDs
entrez_df <- bitr(common_up,
                  fromType = "SYMBOL",
                  toType   = "ENTREZID",
                  OrgDb    = org.Mm.eg.db)
cat("Genes with Entrez ID:", nrow(entrez_df), "\n\n")

# -----------------------------------------------------------------------------
# GO Biological Process enrichment
# -----------------------------------------------------------------------------

go_bp <- enrichGO(
  gene          = entrez_df$ENTREZID,
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

cat("Significant GO BP terms:", nrow(go_bp@result[go_bp@result$p.adjust < 0.05, ]), "\n")

p_go <- dotplot(go_bp, showCategory = 15, font.size = 10) +
  labs(title = "GO Biological Process — genes upregulated in >= 3 organs") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 11))

ggsave("GO_BP_dotplot.pdf", p_go, width = 9, height = 7)
ggsave("GO_BP_dotplot.png", p_go, width = 9, height = 7, dpi = 300)
write.csv(go_bp@result, "GO_BP_results.csv", row.names = FALSE)
cat("Saved: GO_BP_dotplot.pdf, GO_BP_results.csv\n\n")

# -----------------------------------------------------------------------------
# KEGG pathway enrichment
# -----------------------------------------------------------------------------

kegg <- enrichKEGG(
  gene          = entrez_df$ENTREZID,
  organism      = "mmu",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

cat("Significant KEGG pathways:", nrow(kegg@result[kegg@result$p.adjust < 0.05, ]), "\n")

p_kegg <- dotplot(kegg, showCategory = 15, font.size = 10) +
  labs(title = "KEGG pathways — genes upregulated in >= 3 organs") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 11))

ggsave("KEGG_dotplot.pdf", p_kegg, width = 9, height = 7)
ggsave("KEGG_dotplot.png", p_kegg, width = 9, height = 7, dpi = 300)
write.csv(kegg@result, "KEGG_results.csv", row.names = FALSE)
cat("Saved: KEGG_dotplot.pdf, KEGG_results.csv\n")

sessionInfo()
