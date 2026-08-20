# Redundancy Analysis (RDA) of VST-normalized gene expression
#
# Models treatment and organ effects jointly, tests significance with PERMANOVA
# (999 permutations), and produces an ordination plot.
#
# Requires: vst_matrix and metadata in the environment (from 02_counts_to_vst.R)

suppressPackageStartupMessages({
  library(vegan)
  library(ggplot2)
  library(dplyr)
})

treatment_colors <- c(
  "Noinf" = "#7EB07A",
  "Dox"   = "#FF8861",
  "Notra" = "#66A8D6"
)

organ_shapes <- c(
  "Liver"  = 16,
  "Lung"   = 17,
  "Spleen" = 15,
  "Brain"  = 18
)

# -----------------------------------------------------------------------------
# Fit RDA models
# -----------------------------------------------------------------------------

expr_t <- t(vst_matrix)

rda_full      <- rda(expr_t ~ Treatment + Organ, data = metadata)
rda_treatment <- rda(expr_t ~ Treatment,         data = metadata)
rda_organ     <- rda(expr_t ~ Organ,             data = metadata)

var_exp  <- summary(rda_full)$cont$importance[2, 1:2] * 100
r2_treat <- RsquareAdj(rda_treatment)$adj.r.squared
r2_organ <- RsquareAdj(rda_organ)$adj.r.squared

cat("Variance explained by RDA1:", round(var_exp[1], 1), "%\n")
cat("Variance explained by RDA2:", round(var_exp[2], 1), "%\n")
cat("Adj. R² Treatment:", round(r2_treat * 100, 1), "%\n")
cat("Adj. R² Organ:    ", round(r2_organ * 100, 1), "%\n\n")

# -----------------------------------------------------------------------------
# PERMANOVA
# -----------------------------------------------------------------------------

cat("Running PERMANOVA (999 permutations)...\n")

perm_full  <- adonis2(expr_t ~ Treatment * Organ, data = metadata, permutations = 999)
perm_treat <- adonis2(expr_t ~ Treatment,         data = metadata, permutations = 999)
perm_organ <- adonis2(expr_t ~ Organ,             data = metadata, permutations = 999)

cat("\nTreatment: F =", round(perm_treat$F[1], 3),
    "  p =", perm_treat$`Pr(>F)`[1], "\n")
cat("Organ:     F =", round(perm_organ$F[1], 3),
    "  p =", perm_organ$`Pr(>F)`[1], "\n\n")

# -----------------------------------------------------------------------------
# Ordination plot
# -----------------------------------------------------------------------------

scores_df <- data.frame(
  scores(rda_full, display = "sites", scaling = 2),
  Treatment = metadata$Treatment,
  Organ     = metadata$Organ,
  stringsAsFactors = FALSE
)

centroids <- scores_df |>
  group_by(Treatment) |>
  summarise(RDA1 = mean(RDA1), RDA2 = mean(RDA2), .groups = "drop")

p_rda <- ggplot(scores_df, aes(x = RDA1, y = RDA2)) +
  geom_point(aes(color = Treatment, shape = Organ), size = 3.5, alpha = 0.8) +
  geom_point(data = centroids, aes(color = Treatment),
             size = 6, shape = 21, fill = NA, stroke = 2) +
  stat_ellipse(aes(group = Treatment, color = Treatment),
               level = 0.95, linewidth = 1.1) +
  scale_color_manual(values = treatment_colors) +
  scale_shape_manual(values = organ_shapes) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title      = element_text(face = "bold", hjust = 0.5),
    plot.subtitle   = element_text(hjust = 0.5, color = "gray45", size = 10),
    legend.position = "right",
    panel.grid.minor = element_blank()
  ) +
  labs(
    title    = "RDA: VST-normalized gene expression",
    subtitle = paste0(
      "Treatment R² = ", round(r2_treat * 100, 1), "%",
      "  |  Organ R² = ", round(r2_organ * 100, 1), "%"
    ),
    x = paste0("RDA1 (", round(var_exp[1], 1), "%)"),
    y = paste0("RDA2 (", round(var_exp[2], 1), "%)")
  )

ggsave("rda_ordination.pdf", p_rda, width = 8, height = 6)
ggsave("rda_ordination.png", p_rda, width = 8, height = 6, dpi = 300)
cat("Saved: rda_ordination.pdf / .png\n")

write.csv(scores_df, "rda_scores.csv", row.names = FALSE)
cat("Saved: rda_scores.csv\n")

sessionInfo()
