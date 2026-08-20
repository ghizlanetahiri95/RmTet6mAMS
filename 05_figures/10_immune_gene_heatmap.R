# Immune gene expression heatmap (Figure 2I)
#
# Displays Z-score expression of curated immunity-related genes across
# treatment groups (Noinf, Dox, Notra) and organs (Liver, Lung, Spleen, Brain).
# Genes are grouped by functional category.
#
# Requires: vst_symbol and metadata in the environment (from 02_counts_to_vst.R)

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)
})

# -----------------------------------------------------------------------------
# Curated gene panel (matches Figure 2I)
# -----------------------------------------------------------------------------

immune_genes <- list(

  Pro_inflammatory_Cytokines = c(
    "Il1a","Il1b","Il6","Tnf","Il10","Il12a","Il15","Il22",
    "Il23a","Il27","Il17f","Il17a"
  ),
  Chemokines = c(
    "Ccl3","Ccl4","Ccl7","Ccl8","Ccl11","Ccl12","Ccl17","Ccl20","Ccl22",
    "Cxcl1","Cxcl2","Cxcl3","Cxcl5","Cxcl12","Cxcl13"
  ),
  TLRs = c(
    "Tlr1","Tlr2","Tlr3","Tlr4","Tlr5","Tlr6","Tlr8","Tlr9",
    "Tlr11","Tlr12","Tlr13","Ifih1"
  ),
  NLRs_Inflammasome = c(
    "Nlrc3","Nlrc4","Nlrp1a","Nlrp1b","Nlrp3","Nlrp4f",
    "Nlrp6","Nlrp10","Nlrp12","Nod2"
  ),
  CLRs = c(
    "Clec12a","Clec1a","Clec4a2","Clec4b1","Clec4d","Clec4e",
    "Clec4f","Clec4n","Clec5a","Clec7a","Clec10a"
  ),
  Receptors_Cytokine_Signaling = c(
    "Ccr1","Ccr3","Ccr7","Ccrl2","Ifngr1","Ifngr2","Ifnlr1","Tnfr1"
  ),
  Growth_Factors_Mediators = c(
    "Csf2","Csf3","Clcf1","Ptgs2","Hif1a","Egfr","Epha2"
  ),
  Necroptosis_Cell_Death = c(
    "Ripk1","Ripk3","Mlkl","Bach2","Gata3","Lef1","Bcl2"
  )
)

category_colors <- c(
  Pro_inflammatory_Cytokines   = "#E31A1C",
  Chemokines                   = "#FF7F00",
  TLRs                         = "#FDBF6F",
  NLRs_Inflammasome            = "#33A02C",
  CLRs                         = "#1F78B4",
  Receptors_Cytokine_Signaling = "#6A3D9A",
  Growth_Factors_Mediators     = "#B2DF8A",
  Necroptosis_Cell_Death       = "#A6CEE3"
)

# -----------------------------------------------------------------------------
# Build expression matrix: one column per treatment × organ combination
# -----------------------------------------------------------------------------

group_order <- expand.grid(
  Treatment = c("Noinf", "Dox", "Notra"),
  Organ     = c("Liver", "Lung", "Spleen", "Brain"),
  stringsAsFactors = FALSE
)

all_genes <- intersect(unique(unlist(immune_genes)), rownames(vst_symbol))
cat("Genes found in matrix:", length(all_genes), "\n")

expr_group <- sapply(seq_len(nrow(group_order)), function(i) {
  tr <- group_order$Treatment[i]
  org <- group_order$Organ[i]
  idx <- metadata$Treatment == tr & metadata$Organ == org
  if (sum(idx) == 0) return(rep(NA, length(all_genes)))
  rowMeans(vst_symbol[all_genes, idx, drop = FALSE])
})
colnames(expr_group) <- paste(group_order$Treatment, group_order$Organ, sep = "_")
expr_z <- t(scale(t(expr_group)))

# -----------------------------------------------------------------------------
# Annotations
# -----------------------------------------------------------------------------

gene_cat <- sapply(all_genes, function(g) {
  for (cat in names(immune_genes)) if (g %in% immune_genes[[cat]]) return(cat)
  "Other"
})

col_ann <- HeatmapAnnotation(
  Treatment = group_order$Treatment,
  Organ     = group_order$Organ,
  col = list(
    Treatment = c(Noinf = "#7EB07A", Dox = "#FF8861", Notra = "#66A8D6"),
    Organ     = c(Liver = "#D4A017", Lung = "#1F78B4",
                  Spleen = "#6A3D9A", Brain = "#E31A1C")
  ),
  annotation_name_side = "left",
  simple_anno_size = unit(0.4, "cm")
)

row_ann <- rowAnnotation(
  Category = gene_cat,
  col = list(Category = category_colors),
  simple_anno_size = unit(0.35, "cm"),
  annotation_name_side = "top"
)

col_fun <- colorRamp2(c(-2, 0, 2), c("#3093B2", "white", "#CB6677"))

# Order genes by category then by FC Notra vs Noinf in liver
gene_order <- order(
  match(gene_cat, names(immune_genes)),
  -(expr_z[, "Notra_Liver"] - expr_z[, "Noinf_Liver"])
)

ht <- Heatmap(
  expr_z[gene_order, ], name = "Z-score", col = col_fun,
  top_annotation  = col_ann,
  right_annotation = row_ann[gene_order],
  cluster_rows    = FALSE, cluster_columns = FALSE,
  row_split       = gene_cat[gene_order],
  row_title_rot   = 0, row_gap = unit(1.5, "mm"),
  row_names_gp    = gpar(fontsize = 7),
  column_names_gp = gpar(fontsize = 7),
  show_column_names = FALSE,
  border = TRUE
)

pdf("immune_gene_heatmap.pdf", width = 16, height = 14)
draw(ht, heatmap_legend_side = "right")
dev.off()

png("immune_gene_heatmap.png", width = 16, height = 14, units = "in", res = 300)
draw(ht, heatmap_legend_side = "right")
dev.off()

cat("Saved: immune_gene_heatmap.pdf / .png\n")

sessionInfo()
