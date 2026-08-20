# Install all packages required for the analysis

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

bioc_packages <- c(
  "DESeq2",
  "ComplexHeatmap",
  "clusterProfiler",
  "org.Mm.eg.db",
  "AnnotationDbi"
)

cran_packages <- c(
  "ggplot2",
  "vegan",
  "circlize",
  "dplyr",
  "tidyr",
  "ggpubr",
  "scales",
  "gridExtra"
)

cat("Installing CRAN packages...\n")
install.packages(cran_packages, repos = "https://cloud.r-project.org")

cat("Installing Bioconductor packages...\n")
BiocManager::install(bioc_packages, update = FALSE, ask = FALSE)

cat("\nAll packages installed. Checking availability...\n")
all_pkgs <- c(bioc_packages, cran_packages)
for (pkg in all_pkgs) {
  status <- requireNamespace(pkg, quietly = TRUE)
  cat(sprintf("  %-25s %s\n", pkg, ifelse(status, "OK", "FAILED")))
}
