# Record every ggsave() plot object while still writing the requested file.
# Source this before the script under test:
#
#   PLOT_RECORD_DIR=/tmp/plots PLOT_RECORD_BASE=/tmp/results \
#     Rscript -e 'source("analysis/verification/prelude.R"); source("script.R")'
#
# compare_analysis.R renders the recorded objects with a fixed seed.

if (nzchar(Sys.getenv("PLOT_RECORD_DIR"))) {
  .plot_record_dir <- Sys.getenv("PLOT_RECORD_DIR")
  .plot_record_base <- sub("/+$", "", normalizePath(Sys.getenv("PLOT_RECORD_BASE"),
                                                    mustWork = FALSE))
  dir.create(.plot_record_dir, recursive = TRUE, showWarnings = FALSE)

  ggsave <- function(filename, plot = ggplot2::last_plot(), ...,
                     width = NA, height = NA, units = c("in", "cm", "mm"),
                     dpi = 300, limitsize = TRUE, bg = NULL, device = NULL) {
    fn <- normalizePath(filename, mustWork = FALSE)
    rel <- sub(paste0("^", .plot_record_base, "/?"), "", fn)
    rel <- gsub("/+", "/", rel)
    rel <- sub("^/+", "", rel)
    key <- gsub("[^A-Za-z0-9._-]", "_", rel)
    saveRDS(plot, file.path(.plot_record_dir, paste0(key, ".plot.rds")))
    ggplot2::ggsave(
      filename = filename, plot = plot, width = width, height = height,
      units = match.arg(units), dpi = dpi, limitsize = limitsize, bg = bg,
      device = device, ...
    )
  }
}
