# Network discovery: find every res_rc network RDS and annotate its role.

discover_networks <- function(base_dir) {
  cohorts <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)
  cohorts <- cohorts[grepl("^PR", basename(cohorts))]
  rows <- purrr::map(cohorts, function(cohort_dir) {
    files <- list.files(file.path(cohort_dir, "res_rc"),
      pattern = "^g_.+_full_network[.]rds$", full.names = TRUE
    )
    if (!length(files)) return(NULL)
    cohort <- basename(cohort_dir)
    disease <- tail(strsplit(cohort, "_", fixed = TRUE)[[1]], 1)
    condition <- sub("^g_(.+)_full_network[.]rds$", "\\1", basename(files))
    dplyr::tibble(
      cohort = cohort,
      disease_name = disease,
      condition = condition,
      role = dplyr::case_when(
        disease == "ERAS" & condition == "Control" ~ "antibiotic_baseline",
        disease == "ERAS" ~ "antibiotic_followup",
        condition == "Control" ~ "control",
        condition == disease ~ "disease",
        TRUE ~ "descriptive_only"
      ),
      time_order = if (disease == "ERAS") {
        match(condition, c("Control", "D4", "D8", "D42", "D180")) - 1L
      } else {
        rep(NA_integer_, length(condition))
      },
      rds_file = files
    )
  })
  bind_rows_safe(rows) %>% dplyr::arrange(cohort, role, time_order, condition)
}
