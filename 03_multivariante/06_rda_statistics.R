# Permutation tests and post-hoc pairwise PERMANOVA for RDA results
#
# Runs anova.cca on the full RDA model and pairwise PERMANOVA between
# treatment groups. 
#
# Requires: vst_matrix, metadata, and rda_full (from 05_rda_analysis.R)

suppressPackageStartupMessages({
  library(vegan)
  library(dplyr)
})

expr_t <- t(vst_matrix)

# -----------------------------------------------------------------------------
# ANOVA on RDA axes and terms
# -----------------------------------------------------------------------------

cat("ANOVA by axis:\n")
anova_axis  <- anova(rda_full, by = "axis", permutations = 999)
print(anova_axis)

cat("\nANOVA by term:\n")
anova_terms <- anova(rda_full, by = "terms", permutations = 999)
print(anova_terms)

# -----------------------------------------------------------------------------
# Pairwise PERMANOVA between treatment groups
# -----------------------------------------------------------------------------

treatments <- unique(metadata$Treatment)
pairs      <- combn(treatments, 2, simplify = FALSE)

posthoc <- lapply(pairs, function(pair) {
  idx  <- metadata$Treatment %in% pair
  res  <- adonis2(expr_t[idx, ] ~ Treatment,
                  data = metadata[idx, ], permutations = 999)
  data.frame(
    Comparison = paste(pair, collapse = " vs "),
    F_value    = round(res$F[1], 3),
    R2         = round(res$R2[1], 3),
    p_value    = res$`Pr(>F)`[1],
    Significance = dplyr::case_when(
      res$`Pr(>F)`[1] < 0.001 ~ "***",
      res$`Pr(>F)`[1] < 0.01  ~ "**",
      res$`Pr(>F)`[1] < 0.05  ~ "*",
      TRUE                     ~ "ns"
    )
  )
}) |> bind_rows()

cat("\n=== Pairwise PERMANOVA ===\n")
print(posthoc)

write.csv(posthoc, "permanova_posthoc.csv", row.names = FALSE)
cat("Saved: permanova_posthoc.csv\n")

# -----------------------------------------------------------------------------
# Centroid distances
# -----------------------------------------------------------------------------

scores_sites <- scores(rda_full, display = "sites", scaling = 2)
centroids    <- aggregate(scores_sites[, 1:2],
                          by = list(Treatment = metadata$Treatment), mean)
rownames(centroids) <- centroids$Treatment
dist_mat <- as.matrix(dist(centroids[, -1]))

cat("\n=== Centroid distances (RDA space) ===\n")
print(round(dist_mat, 3))

write.csv(as.data.frame(dist_mat), "rda_centroid_distances.csv")
cat("Saved: rda_centroid_distances.csv\n")

sessionInfo()
