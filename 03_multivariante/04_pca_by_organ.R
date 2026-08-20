# PCA of VST-normalized gene expression, stratified by organ
#
# For each organ, runs PCA on the VST matrix, computes PERMANOVA to test
# treatment effects, and exports plots (PDF/PNG/SVG) and scores (CSV).
#
# Requires: vst_matrix and metadata in the environment (from 02_counts_to_vst.R)

suppressPackageStartupMessages({
  library(vegan)
  library(ggplot2)
  library(gridExtra)
  library(dplyr)
})

# Color palette (consistent across all organ panels)
treatment_colors <- c(
  "Noinf" = "#7EB07A",
  "Dox"   = "#FF8861",
  "Notra" = "#66A8D6"
)

# -----------------------------------------------------------------------------
# PCA for a single organ
# -----------------------------------------------------------------------------

pca_single_organ <- function(vst_matrix, metadata, organ) {

  samples <- metadata$Organ == organ
  if (sum(samples) == 0) {
    warning("No samples found for organ: ", organ)
    return(NULL)
  }

  expr   <- t(vst_matrix[, samples])
  meta   <- metadata[samples, ]

  pca_res     <- rda(expr)
  var_exp     <- summary(pca_res)$cont$importance[2, 1:2] * 100
  permanova   <- adonis2(expr ~ Treatment, data = meta, permutations = 999)
  p_val       <- permanova$`Pr(>F)`[1]

  sig_label <- dplyr::case_when(
    p_val < 0.001 ~ "p < 0.001",
    p_val < 0.01  ~ "p < 0.01",
    p_val < 0.05  ~ "p < 0.05",
    TRUE          ~ paste("p =", round(p_val, 3))
  )

  scores_df <- data.frame(
    scores(pca_res, display = "sites", scaling = 2)[, 1:2],
    Treatment = meta$Treatment,
    stringsAsFactors = FALSE
  ) |> setNames(c("PC1", "PC2", "Treatment"))

  centroids <- scores_df |>
    group_by(Treatment) |>
    summarise(PC1 = mean(PC1), PC2 = mean(PC2), .groups = "drop")

  p <- ggplot(scores_df, aes(x = PC1, y = PC2, color = Treatment)) +
    geom_point(size = 3.5, alpha = 0.85) +
    geom_point(data = centroids, size = 5, shape = 21,
               aes(fill = Treatment), color = "white", stroke = 1.5) +
    stat_ellipse(aes(group = Treatment), level = 0.95, linewidth = 1) +
    scale_color_manual(values = treatment_colors) +
    scale_fill_manual(values  = treatment_colors) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title    = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, color = "gray45", size = 10),
      legend.position = "right",
      panel.grid.minor = element_blank()
    ) +
    labs(
      title    = organ,
      subtitle = paste("PERMANOVA:", sig_label),
      x = paste0("PC1 (", round(var_exp[1], 1), "%)"),
      y = paste0("PC2 (", round(var_exp[2], 1), "%)")
    ) +
    guides(fill = "none")

  list(
    organ    = organ,
    plot     = p,
    scores   = scores_df,
    var_exp  = var_exp,
    permanova = permanova
  )
}

# -----------------------------------------------------------------------------
# Run PCA for all organs and assemble 2x2 panel
# -----------------------------------------------------------------------------

organs  <- c("Liver", "Lung", "Spleen", "Brain")
results <- lapply(organs, function(o) pca_single_organ(vst_matrix, metadata, o))
names(results) <- organs

# Print PERMANOVA summary
cat("PERMANOVA results by organ:\n")
for (o in organs) {
  r <- results[[o]]
  if (!is.null(r)) {
    pv <- r$permanova$`Pr(>F)`[1]
    fv <- round(r$permanova$F[1], 3)
    cat(sprintf("  %-8s  F = %.3f  p = %.4f\n", o, fv, pv))
  }
}

# 2x2 panel
panel <- gridExtra::arrangeGrob(
  grobs = lapply(results, `[[`, "plot"),
  ncol = 2
)

ggsave("pca_by_organ.pdf", panel, width = 13, height = 10)
ggsave("pca_by_organ.png", panel, width = 13, height = 10, dpi = 300)
cat("\nSaved: pca_by_organ.pdf / .png\n")

# Export scores per organ
all_scores <- bind_rows(lapply(results, `[[`, "scores"), .id = "Organ")
write.csv(all_scores, "pca_scores_by_organ.csv", row.names = FALSE)
cat("Saved: pca_scores_by_organ.csv\n")

sessionInfo()
