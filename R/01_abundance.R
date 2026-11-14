# Abundance classification and matrix input.

assign_abundance_category <- function(abundance_vector,
                                      high_threshold = 1e-3,
                                      low_threshold_min = 1e-4,
                                      low_threshold_max = 1e-3,
                                      very_low_min = 1e-5,
                                      very_low_max = 1e-4) {
  positive_abundance <- abundance_vector[abundance_vector > 0]

  if (length(positive_abundance) == 0) {
    return("absent")
  }
  categories <- sapply(positive_abundance, function(x) {
    if (x > high_threshold) {
      "high"
    } else if (x > low_threshold_min && x <= low_threshold_max) {
      "low"
    } else if (x > very_low_min && x <= very_low_max) {
      "very_low"
    } else {
      "other"
    }
  })

  freq_table <- table(categories)
  names(freq_table)[which.max(freq_table)]
}

# Collapses the four abundance classes used by the network figures.
abundance_type_from_category <- function(category) {
  ifelse(category == "high", "Rich",
         ifelse(category %in% c("low", "very_low"), "Rare", "Other"))
}

# Keeps the species whose prevalence is at least `threshold`.
# A threshold of 0 disables the filter.
filter_by_prevalence <- function(mat, threshold, present = 0) {
  if (threshold <= 0) {
    return(list(matrix = mat, keep = colnames(mat)))
  }
  prevalence <- colSums(mat > present) / nrow(mat)
  keep <- names(prevalence[prevalence >= threshold])
  list(matrix = mat[, keep, drop = FALSE], keep = keep)
}

merge_and_filter_abundance <- function(file_list, prevalence_threshold = 0.10) {
  cat("Reading", length(file_list), "files...\n")
  matrices <- list()
  for (i in seq_along(file_list)) {
    file <- file_list[i]
    if (file.exists(file)) {
      df <- read.table(file, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
      matrices[[i]] <- as.matrix(df)
      cat("  [", i, "] ", basename(file), ": ", nrow(df), " samples x ",
          ncol(df), " species\n", sep = "")
    } else {
      cat("  Warning: File not found -", file, "\n")
    }
  }

  all_species <- unique(unlist(lapply(matrices, colnames)))
  cat("Total unique species:", length(all_species), "\n")

  merged <- matrix(0,
                   nrow = sum(sapply(matrices, nrow)),
                   ncol = length(all_species),
                   dimnames = list(NULL, all_species))

  row_start <- 1
  for (mat in matrices) {
    n_rows <- nrow(mat)
    rows <- row_start:(row_start + n_rows - 1)
    common_species <- intersect(colnames(mat), all_species)
    merged[rows, common_species] <- mat[, common_species]
    row_start <- row_start + n_rows
  }

  rownames(merged) <- unlist(lapply(matrices, rownames))
  cat("Merged matrix:", nrow(merged), "samples x", ncol(merged), "species\n")

  n_samples <- nrow(merged)
  prevalence <- colSums(merged > 0) / n_samples
  keep_species <- names(prevalence[prevalence >= prevalence_threshold])

  cat("Filtered by", prevalence_threshold * 100, "% prevalence:\n")
  cat("  Before:", ncol(merged), "species\n")
  cat("  After:", length(keep_species), "species\n")

  merged[, keep_species, drop = FALSE]
}
