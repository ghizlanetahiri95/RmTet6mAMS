# Build sample metadata from featureCounts file names


library(dplyr)

# -----------------------------------------------------------------------------
# Option 1: build metadata automatically from file names
# -----------------------------------------------------------------------------

create_metadata_from_files <- function(files_pattern = "*_Counts.txt") {

  files <- list.files(pattern = files_pattern, full.names = FALSE)
  if (length(files) == 0) stop("No files found matching pattern: ", files_pattern)

  cat("Files found:\n")
  print(files)

  sample_names <- gsub("_Counts\\.txt$", "", files)

  parts <- strsplit(sample_names, "_")
  treatment <- sapply(parts, `[`, 1)
  organ_rep  <- sapply(parts, function(x) paste(x[-1], collapse = "_"))
  replicate  <- as.numeric(gsub(".*([0-9]+)$", "\\1", organ_rep))
  organ      <- gsub("[0-9]+$", "", organ_rep)
  organ      <- gsub("_$", "", organ)

  metadata <- data.frame(
    row.names = sample_names,
    Treatment = treatment,
    Organ     = organ,
    Replicate = replicate,
    stringsAsFactors = FALSE
  )

  write.csv(metadata, "metadata.csv", quote = FALSE)
  cat("\nmetadata.csv saved (", nrow(metadata), "samples )\n")
  print(metadata)
  return(metadata)
}

# -----------------------------------------------------------------------------
# Option 2: build metadata manually (use if file names do not follow convention)
# -----------------------------------------------------------------------------

create_metadata_manual <- function() {

  samples <- c(
    "Dox_Liver1",  "Dox_Liver2",  "Dox_Liver3",
    "Dox_Lung1",   "Dox_Lung2",   "Dox_Lung3",
    "Dox_Spleen1", "Dox_Spleen2", "Dox_Spleen3",
    "Dox_Brain1",  "Dox_Brain2",  "Dox_Brain3",
    "Noinf_Liver1",  "Noinf_Liver2",  "Noinf_Liver3",
    "Noinf_Lung1",   "Noinf_Lung2",   "Noinf_Lung3",
    "Noinf_Spleen1", "Noinf_Spleen2", "Noinf_Spleen3",
    "Noinf_Brain1",  "Noinf_Brain2",  "Noinf_Brain3",
    "Notra_Liver1",  "Notra_Liver2",  "Notra_Liver3",
    "Notra_Lung1",   "Notra_Lung2",   "Notra_Lung3",
    "Notra_Spleen1", "Notra_Spleen2", "Notra_Spleen3",
    "Notra_Brain1",  "Notra_Brain2",  "Notra_Brain3"
  )

  metadata <- data.frame(
    row.names = samples,
    Treatment = rep(c("Dox", "Noinf", "Notra"), each = 12),
    Organ     = rep(rep(c("Liver", "Lung", "Spleen", "Brain"), each = 3), 3),
    Replicate = rep(1:3, 12),
    stringsAsFactors = FALSE
  )

  write.csv(metadata, "metadata.csv", quote = FALSE)
  cat("metadata.csv saved (", nrow(metadata), "samples )\n")
  print(metadata)
  return(metadata)
}

# Run: choose one
# metadata <- create_metadata_from_files()
# metadata <- create_metadata_manual()
