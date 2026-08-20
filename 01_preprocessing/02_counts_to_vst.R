# Assemble featureCounts output files into a count matrix and apply
# variance-stabilizing transformation (VST) using DESeq2.
#
# Input:  *_Counts.txt files in the working directory + metadata.csv
# Output: vst_matrix (in environment) + vst_matrix.csv

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
})

# -----------------------------------------------------------------------------
# Parameters
# -----------------------------------------------------------------------------

files_pattern <- "*_Counts.txt"
metadata_file <- "metadata.csv"
min_counts    <- 10          # genes with fewer total counts are excluded

# -----------------------------------------------------------------------------
# Load metadata
# -----------------------------------------------------------------------------

metadata <- read.csv(metadata_file, row.names = 1, stringsAsFactors = FALSE)
cat("Samples:", nrow(metadata), "\n")
cat("Treatments:", paste(unique(metadata$Treatment), collapse = ", "), "\n")
cat("Organs:",    paste(unique(metadata$Organ),     collapse = ", "), "\n\n")

# -----------------------------------------------------------------------------
# Read and merge featureCounts files
# -----------------------------------------------------------------------------

read_counts <- function(path) {
  df <- read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  sample_name <- gsub("_Counts\\.txt$", "", basename(path))
  data.frame(Geneid = df[, 1], counts = df[, 2],
             row.names = NULL, stringsAsFactors = FALSE) |>
    setNames(c("Geneid", sample_name))
}

files <- list.files(pattern = files_pattern, full.names = TRUE)
if (length(files) == 0) stop("No count files found.")

counts_list   <- lapply(files, read_counts)
counts_merged <- Reduce(function(a, b) merge(a, b, by = "Geneid", all = TRUE),
                        counts_list)

rownames(counts_merged) <- counts_merged$Geneid
counts_raw <- as.matrix(counts_merged[, -1])
counts_raw[is.na(counts_raw)] <- 0

cat("Count matrix dimensions:", nrow(counts_raw), "genes x", ncol(counts_raw), "samples\n")
cat("Mean library size:", format(round(mean(colSums(counts_raw))), big.mark = ","), "reads\n\n")

# -----------------------------------------------------------------------------
# DESeq2: filter low-count genes and apply VST
# -----------------------------------------------------------------------------

# Align sample order between counts and metadata
metadata <- metadata[colnames(counts_raw), ]

dds <- DESeqDataSetFromMatrix(
  countData = counts_raw,
  colData   = metadata,
  design    = ~ Treatment + Organ
)

keep <- rowSums(counts(dds)) >= min_counts
dds  <- dds[keep, ]
cat("Genes retained after filtering:", nrow(dds), "\n")

vst_data   <- vst(dds, blind = TRUE)
vst_matrix <- assay(vst_data)

cat("VST range:", paste(round(range(vst_matrix), 2), collapse = " to "), "\n\n")

# -----------------------------------------------------------------------------
# Convert ENSEMBL IDs to gene symbols
# -----------------------------------------------------------------------------

library(org.Mm.eg.db)
library(AnnotationDbi)

ensembl_ids <- rownames(vst_matrix)
gene_symbols <- mapIds(org.Mm.eg.db,
                       keys      = ensembl_ids,
                       column    = "SYMBOL",
                       keytype   = "ENSEMBL",
                       multiVals = "first")

keep_sym   <- !is.na(gene_symbols) & !duplicated(gene_symbols)
vst_symbol <- vst_matrix[keep_sym, ]
rownames(vst_symbol) <- gene_symbols[keep_sym]

cat("Genes with SYMBOL annotation:", nrow(vst_symbol), "\n\n")

# -----------------------------------------------------------------------------
# Save
# -----------------------------------------------------------------------------

write.csv(vst_matrix, "vst_matrix_ensembl.csv")
write.csv(vst_symbol, "vst_matrix_symbol.csv")
cat("Saved: vst_matrix_ensembl.csv and vst_matrix_symbol.csv\n")

sessionInfo()
