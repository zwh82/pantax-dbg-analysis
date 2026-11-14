# Configuration helpers: cohort table, output paths, and validation.

read_cohort_table <- function(path, merged = FALSE) {
  config_data <- read.table(path, sep = "\t", header = FALSE,
                            stringsAsFactors = FALSE, comment.char = "#")
  if (nrow(config_data) == 0) {
    stop("Cohort table is empty: ", path)
  }
  if (merged) {
    colnames(config_data) <- c("disease", "metadata_file", "species_abund",
                               "maaslin3_res_file", "project_name",
                               "all_read_count", "all_species_abund")
  } else if (ncol(config_data) == 6) {
    colnames(config_data) <- c("disease", "cohort", "metadata_file",
                               "species_relative_abund", "species_abund",
                               "maaslin3_res_file")
  } else if (ncol(config_data) == 5) {
    colnames(config_data) <- c("disease", "cohort", "metadata_file",
                               "species_relative_abund", "maaslin3_res_file")
  } else {
    stop("Unsupported cohort table layout (", ncol(config_data),
         " columns): ", path)
  }
  config_data
}

resolve_output_dir <- function(cfg) {
  base <- cfg$output$dir
  if (!grepl("/$", base)) base <- paste0(base, "/")
  if (cfg$runs$network_only) {
    return(paste0(base, "all_country/"))
  }
  if (cfg$runs$core_network_only) {
    return(paste0(base, "all_country_core/"))
  }
  base
}

resolve_output_suffix <- function(cfg) {
  if (!is.null(cfg$output$suffix)) {
    return(cfg$output$suffix)
  }

  read_count <- identical(cfg$input$data_type, "read_count")
  if (cfg$network$filter_by_weight) {
    suffix <- if (read_count) "res_rc_filter" else "res_raclr_filter"
  } else if (cfg$network$filter_by_confidence) {
    suffix <- "res_raclr_econ_filter"
  } else if (read_count) {
    suffix <- if (cfg$input$fill_zero_counts) "res_rc_fill" else "res_rc"
  } else {
    suffix <- "res_raclr_no_filter"
  }

  if (identical(cfg$network$tool, "flashweave")) {
    suffix <- if (read_count) {
      "res_flashweave_rc"
    } else {
      paste0("res_flashweave_", cfg$input$prevalence)
    }
  } else if (identical(cfg$network$tool, "sparcc")) {
    suffix <- if (read_count) "res_sparcc_rc" else "res_sparcc"
  }

  if (cfg$input$floor > 0) {
    suffix <- paste0(suffix, "_FilterLow-", cfg$input$floor)
  }
  suffix
}

validate_config <- function(cfg) {
  required <- c("input", "network", "abundance", "hub", "output", "runs", "runtime")
  missing <- setdiff(required, names(cfg))
  if (length(missing) > 0) {
    stop("Config is missing sections: ", paste(missing, collapse = ", "))
  }
  if (!cfg$input$data_type %in% c("read_count", "relative_abundance")) {
    stop("input$data_type must be 'read_count' or 'relative_abundance'")
  }
  if (!cfg$network$tool %in% c("spiec-easi", "flashweave", "sparcc")) {
    stop("network$tool must be 'spiec-easi', 'flashweave', or 'sparcc'")
  }
  invisible(TRUE)
}
