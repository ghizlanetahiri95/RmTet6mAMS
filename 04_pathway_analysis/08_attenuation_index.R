# Doxycycline Attenuation Index for RCD, Inflammation, and Migration
#
# Quantifies the degree to which doxycycline treatment reverses
# infection-induced transcriptional changes for curated gene sets.
#
# Attenuation Index (AI) = 1 - (FC_Dox/Noinf / FC_Notra/Noinf)
#   AI = 1: complete reversal to non-infected baseline
#   AI = 0: no attenuation
#
# Requires: vst_symbol and metadata in the environment (from 02_counts_to_vst.R)

suppressPackageStartupMessages({
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)
})

# -----------------------------------------------------------------------------
# Gene sets
# -----------------------------------------------------------------------------

gene_sets <- list(

  RCD = list(
    Necroptosis = c("Ripk1", "Ripk3", "Mlkl", "Tnfrsf1a", "Zbp1"),
    Pyroptosis  = c("Gsdmd", "Casp1", "Casp4", "Casp11", "Nlrp3",
                    "Nlrc4", "Aim2", "Pycard", "Il1b", "Il18"),
    Apoptosis   = c("Casp3", "Casp8", "Casp9", "Bax", "Bid", "Fas", "Fadd")
  ),

  Inflammation = list(
    Cytokines      = c("Il1a", "Il1b", "Il6", "Tnf", "Il12a", "Il12b", "Il18", "Il23a"),
    NFkB_signaling = c("Nfkb1", "Nfkb2", "Rela", "Relb", "Nfkbia", "Ikbkb", "Ikbkg"),
    Inflammasome   = c("Nlrp3", "Nlrc4", "Pycard", "Casp1", "Il1b", "Il18"),
    Acute_phase    = c("Saa1", "Saa2", "Saa3", "Hp", "Crp", "Orm1")
  ),

  Migration = list(
    Chemokines  = c("Cxcl1", "Cxcl2", "Cxcl3", "Cxcl5", "Ccl2", "Ccl3", "Ccl4", "Ccl7"),
    Adhesion    = c("Icam1", "Vcam1", "Sele", "Selp", "Sell", "Itgam", "Itgax", "Itgb2"),
    Recruitment = c("Csf3", "Csf2", "S100a8", "S100a9", "Lbp")
  )
)

# -----------------------------------------------------------------------------
# Calculate Attenuation Index
# -----------------------------------------------------------------------------

calc_ai <- function(vst, meta, gene_set, organ = NULL) {

  if (!is.null(organ)) {
    idx <- meta$Organ == organ
    vst  <- vst[,  idx]
    meta <- meta[idx, ]
  }

  s_noinf <- meta$Treatment == "Noinf"
  s_dox   <- meta$Treatment == "Dox"
  s_notra <- meta$Treatment == "Notra"

  results <- lapply(names(gene_set), function(main) {
    lapply(names(gene_set[[main]]), function(sub) {
      genes <- intersect(gene_set[[main]][[sub]], rownames(vst))
      if (length(genes) == 0) return(NULL)

      e_noinf <- mean(colMeans(vst[genes, s_noinf, drop = FALSE]))
      e_dox   <- mean(colMeans(vst[genes, s_dox,   drop = FALSE]))
      e_notra <- mean(colMeans(vst[genes, s_notra, drop = FALSE]))

      fc_notra <- e_notra - e_noinf
      fc_dox   <- e_dox   - e_noinf

      ai <- if (fc_notra != 0) max(0, min(1, 1 - fc_dox / fc_notra)) else 0

      data.frame(
        Process = main, Sub_process = sub,
        N_genes = length(genes),
        FC_Notra = fc_notra, FC_Dox = fc_dox,
        AI = ai,
        stringsAsFactors = FALSE
      )
    }) |> bind_rows()
  }) |> bind_rows()

  results
}

# Overall attenuation
ai_overall <- calc_ai(vst_symbol, metadata, gene_sets)
cat("=== Overall Attenuation Index ===\n")
print(ai_overall[order(-ai_overall$AI), c("Process","Sub_process","AI","FC_Notra","FC_Dox")])
write.csv(ai_overall, "attenuation_index_overall.csv", row.names = FALSE)

# -----------------------------------------------------------------------------
# Barplot: overall AI by sub-process
# -----------------------------------------------------------------------------

process_colors <- c(RCD = "#E31A1C", Inflammation = "#FF7F00", Migration = "#377EB8")

ai_plot <- ai_overall |>
  mutate(label = paste(Process, Sub_process, sep = "\n")) |>
  arrange(desc(AI)) |>
  mutate(label = factor(label, levels = label))

p_bar <- ggplot(ai_plot, aes(x = label, y = AI * 100, fill = Process)) +
  geom_col(alpha = 0.85) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "gray40") +
  scale_fill_manual(values = process_colors) +
  scale_y_continuous(limits = c(0, 100)) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        panel.grid.major.x = element_blank(),
        legend.position = "bottom") +
  labs(x = NULL, y = "Attenuation Index (%)", fill = "Process",
       title = "Doxycycline attenuation of infection-induced transcriptional responses")

ggsave("attenuation_barplot.pdf", p_bar, width = 10, height = 6)
ggsave("attenuation_barplot.png", p_bar, width = 10, height = 6, dpi = 300)
cat("Saved: attenuation_barplot.pdf\n")

# -----------------------------------------------------------------------------
# Heatmap: treatment gradient (Noinf → Dox → Notra)
# -----------------------------------------------------------------------------

# Select top genes per sub-process (ranked by FC in Notra vs Noinf)
top_n <- 5
selected <- lapply(names(gene_sets), function(main) {
  lapply(names(gene_sets[[main]]), function(sub) {
    genes <- intersect(gene_sets[[main]][[sub]], rownames(vst_symbol))
    if (length(genes) == 0) return(NULL)
    s_noinf <- metadata$Treatment == "Noinf"
    s_notra <- metadata$Treatment == "Notra"
    fc <- rowMeans(vst_symbol[genes, s_notra, drop=FALSE]) -
          rowMeans(vst_symbol[genes, s_noinf, drop=FALSE])
    top <- names(sort(fc, decreasing = TRUE))[seq_len(min(top_n, length(fc)))]
    setNames(list(top), paste(main, sub, sep = "_"))
  })
}) |> unlist(recursive = FALSE)
selected <- Filter(Negate(is.null), selected)

all_genes <- unique(unlist(selected))

treatments <- c("Noinf", "Dox", "Notra")
expr_avg <- sapply(treatments, function(tr)
  rowMeans(vst_symbol[all_genes, metadata$Treatment == tr, drop = FALSE]))
expr_z <- t(scale(t(expr_avg)))

gene_cat <- sapply(all_genes, function(g) {
  for (n in names(selected)) if (g %in% selected[[n]]) return(strsplit(n, "_")[[1]][1])
})

col_ann <- HeatmapAnnotation(
  Treatment = factor(treatments, levels = treatments),
  col = list(Treatment = c(Noinf = "#7EB07A", Dox = "#FF8861", Notra = "#66A8D6")),
  annotation_name_side = "left", simple_anno_size = unit(0.5, "cm")
)
row_ann <- rowAnnotation(
  Process = gene_cat,
  col = list(Process = process_colors),
  simple_anno_size = unit(0.35, "cm"), annotation_name_side = "top"
)

col_fun <- colorRamp2(c(-2, 0, 2), c("#3093B2", "white", "#CB6677"))

o <- order(match(gene_cat, names(process_colors)),
           -(expr_z[, "Notra"] - expr_z[, "Noinf"]))

ht <- Heatmap(
  expr_z[o, ], name = "Z-score", col = col_fun,
  top_annotation = col_ann,
  right_annotation = row_ann[o],
  cluster_rows = FALSE, cluster_columns = FALSE,
  row_split = gene_cat[o], row_title_rot = 0, row_gap = unit(2, "mm"),
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 11, fontface = "bold"),
  column_title = "Treatment response gradient",
  column_title_gp = gpar(fontsize = 11, fontface = "bold"),
  border = TRUE, width = unit(5, "cm")
)

pdf("attenuation_heatmap.pdf", width = 8, height = 13)
draw(ht, heatmap_legend_side = "right")
dev.off()
cat("Saved: attenuation_heatmap.pdf\n")

# -----------------------------------------------------------------------------
# Barplot: organ-specific attenuation
# -----------------------------------------------------------------------------

ai_organ <- bind_rows(lapply(c("Liver","Lung","Spleen","Brain"), function(org) {
  bind_rows(lapply(names(gene_sets), function(main) {
    genes <- intersect(unique(unlist(gene_sets[[main]])), rownames(vst_symbol))
    if (length(genes) == 0) return(NULL)
    r <- calc_ai(vst_symbol, metadata, setNames(list(setNames(list(genes), main)), main), organ = org)
    r$Organ <- org
    r
  }))
}))

p_organ <- ggplot(ai_organ, aes(x = Organ, y = AI * 100, fill = Process)) +
  geom_col(position = "dodge", alpha = 0.85) +
  scale_fill_manual(values = process_colors) +
  scale_y_continuous(limits = c(0, 100)) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.x = element_blank(), legend.position = "bottom") +
  labs(x = NULL, y = "Attenuation Index (%)", fill = "Process",
       title = "Organ-specific doxycycline attenuation")

ggsave("attenuation_by_organ.pdf", p_organ, width = 8, height = 6)
ggsave("attenuation_by_organ.png", p_organ, width = 8, height = 6, dpi = 300)
cat("Saved: attenuation_by_organ.pdf\n")

sessionInfo()
