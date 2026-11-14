#!/usr/bin/env Rscript
# Rebuild the Hadza metadata tables for the 179-sample matrices.
#
# Usage:
#   Rscript config/build_hadza_metadata.R \
#     /snode_data/jialu/hazha/Hadza_Seasonality_Expansion_179_20261113/tables/sample_manifest.tsv \
#     /snode_data/jialu/hazha/species/PRJEB49206_278_hazha_species_abundance_matrix.txt \
#     config
#
# Fine-season labels follow the collection date, matching the original
# five-run metadata:
#   Jan 2014 -> 2014-EW (Control); Feb 2014 -> 2014-LW
#   May-Jul 2014 -> 2014-ED; Aug 2014 -> 2014-LD
#   Sep-Oct 2013 -> 2013-LD
# The three-timepoint table uses Season directly.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: build_hadza_metadata.R MANIFEST MATRIX OUTDIR")
}
manifest_file <- args[1]
matrix_file <- args[2]
outdir <- args[3]

read_tsv_chr <- function(path) {
  read.delim(path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE,
    colClasses = "character", na.strings = character(0))
}

manifest <- read_tsv_chr(manifest_file)
matrix_ids <- data.table::fread(matrix_file, select = 1L, header = TRUE,
  data.table = FALSE)[[1]]
matrix_ids <- trimws(as.character(matrix_ids))

missing <- setdiff(matrix_ids, manifest$Sample_ID)
if (length(missing)) {
  stop("No manifest row for matrix samples: ", paste(missing, collapse = ", "))
}
# Keep the matrix order so the metadata and matrices line up one to one.
m <- manifest[match(matrix_ids, manifest$Sample_ID), , drop = FALSE]

empty_to_na <- function(x) {
  x <- trimws(as.character(x))
  x[!nzchar(x)] <- "NA"
  x
}

format_age <- function(x) {
  x <- trimws(as.character(x))
  out <- ifelse(x %in% c("", "NA"), "NA", sub("\\.0+$", "", x))
  out
}

format_date <- function(x) {
  d <- as.Date(x)
  if (any(is.na(d))) stop("Unparseable collection_date in manifest.")
  paste(format(d, "%Y"), as.integer(format(d, "%m")), as.integer(format(d, "%d")),
    sep = "/")
}

sample_date <- as.Date(m$collection_date)
year <- as.integer(format(sample_date, "%Y"))
month <- as.integer(format(sample_date, "%m"))
fine_season <- ifelse(year == 2013 & month %in% 9:10, "2013-LD",
  ifelse(year == 2014 & month == 1, "2014-EW",
    ifelse(year == 2014 & month == 2, "2014-LW",
      ifelse(year == 2014 & month %in% 5:7, "2014-ED",
        ifelse(year == 2014 & month == 8, "2014-LD", NA_character_)
      )
    )
  )
)
if (any(is.na(fine_season))) {
  stop("Unmapped collection date: ",
    paste(unique(m$collection_date[is.na(fine_season)]), collapse = "; "))
}

age <- format_age(m$age)
age_num <- suppressWarnings(as.numeric(age))
adult <- ifelse(is.na(age_num), "False", ifelse(age_num >= 18, "True", "False"))

base <- data.frame(
  Sample_ID = m$Sample_ID,
  Sample_Name = m$library_id,
  Population = m$author_population,
  Season = m$original_season_label,
  Subject_ID = m$donor_id,
  Collection_date = format_date(m$collection_date),
  Reads_accession = m$Run_accession,
  BioSample = m$BioSample,
  age_years = age,
  camp = empty_to_na(m$camp),
  fine_season = fine_season,
  profile_complete = "True",
  known_adult_ge18 = adult,
  stringsAsFactors = FALSE
)

five <- base
five$Cohort <- "Hadza_five_timepoints"
five$Standard_Disease <- ifelse(five$fine_season == "2014-EW", "Control",
  five$fine_season)

three <- base
three$Cohort <- "Hadza_three_timepoints"
three$Standard_Disease <- ifelse(three$Season == "Wet2014", "Control",
  three$Season)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
write.table(five, file.path(outdir, "hadza_metadata_five.tsv"), sep = "\t",
  quote = FALSE, row.names = FALSE, na = "NA")
write.table(three, file.path(outdir, "hadza_metadata_three.tsv"), sep = "\t",
  quote = FALSE, row.names = FALSE, na = "NA")

cat("Samples:", nrow(five), "\n")
cat("Five-timepoint groups:\n")
print(table(five$Standard_Disease))
cat("Three-timepoint groups:\n")
print(table(three$Standard_Disease))
