# Shared helpers for the analysis / plotting scripts.

ensure_dir <- function(path) {
  dir.create(path, showWarnings = FALSE, recursive = TRUE)
  path
}

resolve_result_path <- function(cfg, path) {
  if (is.null(path) || !nzchar(path)) return(NULL)
  if (grepl("^/", path)) path else file.path(cfg$results_dir, path)
}

# Cohort networks and analysis outputs can live in separate roots.  When
# `network_dir` is unset the legacy layout (both under results_dir) is used.
network_dir <- function(cfg) {
  if (!is.null(cfg$network_dir) && nzchar(cfg$network_dir)) {
    cfg$network_dir
  } else {
    cfg$results_dir
  }
}

analysis_suffixes <- function(cfg) {
  suffixes <- cfg$suffix
  if (is.null(suffixes) || length(suffixes) == 0) {
    stop("analysis config: `suffix` must name at least one result directory")
  }
  as.character(suffixes)
}

validate_analysis_config <- function(cfg) {
  if (!dir.exists(cfg$results_dir)) {
    stop("results_dir does not exist: ", cfg$results_dir)
  }
  if (!dir.exists(network_dir(cfg))) {
    stop("network_dir does not exist: ", network_dir(cfg))
  }
  analysis_suffixes(cfg)
  if (!is.null(cfg$disease_map_file) && !file.exists(cfg$disease_map_file)) {
    stop("disease_map_file does not exist: ", cfg$disease_map_file)
  }
  invisible(TRUE)
}

# Cohorts are the immediate sub-directories whose name starts with "PR".
# The disease label is the last "_"-separated token of the directory name.
get_disease_cohorts <- function(base_dir) {
  all_dirs <- list.dirs(base_dir, recursive = FALSE, full.names = FALSE)
  all_dirs <- all_dirs[grepl("^PR", all_dirs)]
  diseases <- vapply(strsplit(all_dirs, "_"), function(x) tail(x, 1), character(1))
  split(all_dirs, diseases)
}

read_disease_map <- function(cfg) {
  path <- cfg$disease_map_file
  if (is.null(path)) {
    path <- file.path(cfg$results_dir, "scripts", "disease_name_mapping.tsv")
  }
  if (!file.exists(path)) {
    stop("disease name mapping not found: ", path)
  }
  read.delim(
    path,
    header = FALSE,
    col.names = c("Disease_old", "Disease_new"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

replace_disease_name <- function(df, disease_map) {
  df %>%
    left_join(disease_map, by = c("Disease" = "Disease_old")) %>%
    mutate(Disease = ifelse(is.na(Disease_new), Disease, Disease_new)) %>%
    select(-Disease_new)
}

replace_disease_name2 <- function(df, disease_map, cols = c("Disease")) {
  for (col in cols) {
    if (col %in% colnames(df)) {
      old_value <- as.character(df[[col]])
      new_value <- disease_map$Disease_new[match(old_value, disease_map$Disease_old)]
      df[[col]] <- ifelse(!is.na(new_value), new_value, old_value)
    }
  }
  df
}

replace_single_disease_name <- function(disease, disease_map) {
  disease_new <- disease_map$Disease_new[match(disease, disease_map$Disease_old)]
  ifelse(is.na(disease_new), disease, disease_new)
}

# The disease label used to find a cohort's disease network file.
cohort_disease_token <- function(cohort) {
  tail(strsplit(cohort, "_")[[1]], 1)
}

# First row of a one-row data frame or the row whose Group matches `disease`.
pick_disease_row <- function(rows, disease, cohort) {
  if (nrow(rows) == 1) return(rows)
  exact <- rows %>% filter(toupper(Group) == toupper(disease))
  if (nrow(exact) > 0) return(exact)
  warning("No ", disease, " for cohort ", cohort, " - using ", rows$Group[1])
  rows[1, , drop = FALSE]
}

# Keep only the y-axis label of the first panel in each row of a patchwork
# grid.  Used by the multi-disease rho figures.
keep_first_title_each_row <- function(plot_list, n_col) {
  lapply(seq_along(plot_list), function(i) {
    if ((i - 1) %% n_col == 0) {
      plot_list[[i]]
    } else {
      plot_list[[i]] + labs(y = NULL)
    }
  })
}

# Replace "Rich" with "Abundant" in edge labels.
rename_rich_to_abundant <- function(x) {
  str_replace_all(as.character(x), "Rich", "Abundant")
}
