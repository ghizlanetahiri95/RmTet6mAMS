# GM-CSF signaling analysis
#
# Tests whether infection suppresses GM-CSF signaling and whether
# doxycycline treatment restores it. Genes are classified by expression
# pattern across the three treatment groups (Noinf, Dox, Notra).
# Pathway-level evidence is provided by GSEA on GO Biological Process terms.
#
# Requires: vst_symbol and metadata in the environment (from 02_counts_to_vst.R)

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(ggplot2)
  library(dplyr)
})

# -----------------------------------------------------------------------------
# Curated GM-CSF gene sets
# -----------------------------------------------------------------------------

gmcsf_genes <- list(
  Production = c("Csf2"),
  Receptor   = c("Csf2ra", "Csf2rb", "Csf2rb2"),
  JAK_STAT   = c("Jak1", "Jak2", "Stat5a", "Stat5b", "Stat3"),
  Myeloid_TF = c("Cebpa", "Cebpb", "Spi1", "Irf8", "Irf4"),
  Survival   = c("Mcl1", "Bcl2", "Bcl2l1"),
  Phagocyte  = c("Fcgr1", "Fcgr2b", "Fcgr3", "Itgam", "Itgax",
                 "Cybb", "Ncf1", "Ncf2")
)

pathway_colors <- c(
  Production = "#E31A1C", Receptor  = "#FF7F00", JAK_STAT   = "#6A3D9A",
  Myeloid_TF = "#33A02C", Survival  = "#1F78B4", Phagocyte  = "#B15928"
)

# -----------------------------------------------------------------------------
# Expression averages and gene-level summary
# -----------------------------------------------------------------------------

genes_found <- intersect(unique(unlist(gmcsf_genes)), rownames(vst_symbol))
cat("GM-CSF genes found:", length(genes_found), "of", length(unique(unlist(gmcsf_genes))), "\n")

treatments <- c("Noinf", "Dox", "Notra")
expr_avg <- sapply(treatments, function(tr)
  rowMeans(vst_symbol[genes_found, metadata$Treatment == tr, drop = FALSE]))
expr_z <- t(scale(t(expr_avg)))

gene_cat <- sapply(genes_found, function(g) {
  for (cat in names(gmcsf_genes)) if (g %in% gmcsf_genes[[cat]]) return(cat)
})

summary_df <- data.frame(
  Gene     = genes_found,
  Category = gene_cat,
  FC_Notra = expr_avg[genes_found, "Notra"] - expr_avg[genes_found, "Noinf"],
  FC_Dox   = expr_avg[genes_found, "Dox"]   - expr_avg[genes_found, "Noinf"],
  stringsAsFactors = FALSE
) |>
  mutate(Pattern = case_when(
    FC_Notra < -0.3 & FC_Dox > FC_Notra ~ "Suppressed_Restored",
    FC_Notra < -0.3                      ~ "Suppressed_NotRestored",
    FC_Notra >  0.3                      ~ "Induced",
    TRUE                                 ~ "NoChange"
  ))

cat("\n=== Gene-level summary (sorted by infection FC) ===\n")
print(summary_df[order(summary_df$FC_Notra),
                 c("Gene","Category","FC_Notra","FC_Dox","Pattern")])

write.csv(summary_df, "GMCSF_gene_summary.csv", row.names = FALSE)

# -----------------------------------------------------------------------------
# Heatmap
# -----------------------------------------------------------------------------

gene_order <- order(match(gene_cat, names(gmcsf_genes)),
                    summary_df$FC_Notra)

col_ann <- HeatmapAnnotation(
  Treatment = factor(treatments, levels = treatments),
  col = list(Treatment = c(Noinf = "#7EB07A", Dox = "#FF8861", Notra = "#66A8D6")),
  annotation_name_side = "left", simple_anno_size = unit(0.5, "cm")
)
row_ann <- rowAnnotation(
  Pathway = gene_cat[gene_order],
  col = list(Pathway = pathway_colors),
  simple_anno_size = unit(0.35, "cm"), annotation_name_side = "top"
)

col_fun <- colorRamp2(c(-2, 0, 2), c("#3093B2", "white", "#CB6677"))

ht <- Heatmap(
  expr_z[gene_order, ], name = "Z-score", col = col_fun,
  top_annotation = col_ann,
  right_annotation = row_ann,
  cluster_rows = FALSE, cluster_columns = FALSE,
  row_split = gene_cat[gene_order], row_title_rot = 0,
  row_gap = unit(2, "mm"),
  row_names_gp = gpar(fontsize = 9),
  column_names_gp = gpar(fontsize = 11, fontface = "bold"),
  column_title = "GM-CSF signaling across treatment conditions",
  column_title_gp = gpar(fontsize = 11, fontface = "bold"),
  border = TRUE, width = unit(5, "cm")
)

pdf("GMCSF_heatmap.pdf", width = 8, height = 10)
draw(ht, heatmap_legend_side = "right")
dev.off()
cat("Saved: GMCSF_heatmap.pdf\n")

# -----------------------------------------------------------------------------
# GSEA — pathways suppressed or activated during infection
# -----------------------------------------------------------------------------

fc_all <- rowMeans(vst_symbol[, metadata$Treatment == "Notra", drop=FALSE]) -
          rowMeans(vst_symbol[, metadata$Treatment == "Noinf", drop=FALSE])
fc_ranked <- sort(fc_all, decreasing = TRUE)

gene_df   <- bitr(names(fc_ranked), fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Mm.eg.db)
fc_entrez <- fc_ranked[gene_df$SYMBOL]
names(fc_entrez) <- gene_df$ENTREZID
fc_entrez <- sort(fc_entrez, decreasing = TRUE)

gsea_res <- gseGO(
  geneList     = fc_entrez,
  OrgDb        = org.Mm.eg.db,
  ont          = "BP",
  minGSSize    = 10,
  maxGSSize    = 500,
  pvalueCutoff = 0.05,
  verbose      = FALSE
)

write.csv(gsea_res@result, "GMCSF_GSEA_full.csv", row.names = FALSE)

# Subset to immune / myeloid / phagocyte terms
keywords <- c("granulocyte","neutrophil","phagocyte","myeloid","CSF",
              "JAK","STAT","monocyte","macrophage","innate","leukocyte")

gsea_immune <- gsea_res@result |>
  filter(grepl(paste(keywords, collapse="|"), Description, ignore.case=TRUE)) |>
  arrange(NES)

cat("\n=== Immune-related GSEA terms (NES < 0 = suppressed in infection) ===\n")
print(gsea_immune[, c("Description","NES","p.adjust")])

# Dotplot of top terms with NES < 0
gsea_neg        <- gsea_res
gsea_neg@result <- gsea_res@result[gsea_res@result$NES < 0, ]
if (nrow(gsea_neg@result) > 0) {
  pdf("GMCSF_GSEA_suppressed.pdf", width = 11, height = 8)
  print(dotplot(gsea_neg, showCategory = 20,
                title = "Pathways suppressed during infection (NES < 0)"))
  dev.off()
  cat("Saved: GMCSF_GSEA_suppressed.pdf\n")
}

sessionInfo()
