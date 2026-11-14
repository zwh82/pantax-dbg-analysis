# Config file loader shared by run_pipeline.R and analysis/run_analysis.R.
#
# YAML files are the default.  R files are still accepted for existing
# recipes: they must define `config` (pipeline) or `analysis_config`
# (analysis).  YAML sequences become atomic vectors; YAML mappings stay
# named lists (use unlist() where a named vector is required).

normalize_config <- function(x) {
  if (!is.list(x)) return(x)
  x <- lapply(x, normalize_config)
  unnamed <- is.null(names(x)) || all(!nzchar(names(x)))
  if (unnamed && length(x) == 0) return(character(0))
  if (unnamed && all(vapply(x, function(e) is.atomic(e) && length(e) == 1,
                            logical(1)))) {
    return(unlist(x, use.names = FALSE))
  }
  x
}

read_config_file <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("yaml", "yml")) {
    if (!requireNamespace("yaml", quietly = TRUE)) {
      stop("Reading a YAML config requires the 'yaml' package; ",
           "use an .R config file or install yaml")
    }
    cfg <- normalize_config(yaml::read_yaml(path))
  } else {
    env <- new.env(parent = baseenv())
    sys.source(path, envir = env)
    cfg <- if (!is.null(env$config)) env$config else env$analysis_config
  }
  if (!is.list(cfg)) stop("Config is empty or malformed: ", path)
  cfg
}

# yaml::as.yaml() writes named atomic vectors as unnamed sequences, so convert
# them to named lists first; the result is a YAML mapping with the same names.
config_for_yaml <- function(x) {
  if (is.list(x)) return(lapply(x, config_for_yaml))
  if (!is.null(names(x)) && any(nzchar(names(x)))) return(as.list(x))
  x
}

write_config_file <- function(cfg, path) {
  ext <- tolower(tools::file_ext(path))
  if (dirname(path) != ".") {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  }
  if (ext %in% c("yaml", "yml")) {
    if (!requireNamespace("yaml", quietly = TRUE)) {
      stop("Writing a YAML config requires the 'yaml' package")
    }
    yaml::write_yaml(config_for_yaml(cfg), path)
  } else {
    con <- file(path, "w")
    writeLines("config <-", con)
    dput(cfg, con)
    close(con)
  }
  invisible(path)
}
