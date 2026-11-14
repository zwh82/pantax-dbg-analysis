# Hub overlap with differential species markers.

analyze_hub_differential_overlap <- function(hub_results, group_name, sig_df,
                                             outdir, outflag = "") {
  cat("\n========== Analyzing Hub-Differential Overlap for", group_name, "==========\n")

  consensus_hubs <- hub_results$consensus_hubs

  if (length(consensus_hubs) == 0) {
    cat("No consensus hubs found in", group_name, "\n")
    return(NULL)
  }

  hub_differential <- sig_df[sig_df$species_taxid %in% consensus_hubs, ]

  n_hubs <- length(consensus_hubs)
  n_diff_in_hubs <- nrow(hub_differential)
  n_risk_in_hubs <- sum(hub_differential$sig_marker == "risk")
  n_protective_in_hubs <- sum(hub_differential$sig_marker == "protective")

  cat("\nHub nodes in", group_name, ":\n")
  cat("  Total hubs:", n_hubs, "\n")
  cat("  Differential markers in hubs:", n_diff_in_hubs,
      sprintf("(%.1f%%)", n_diff_in_hubs / n_hubs * 100), "\n")
  cat("    - Risk hubs:", n_risk_in_hubs, "\n")
  cat("    - Protective hubs:", n_protective_in_hubs, "\n")

  result_df <- data.frame(
    Group = group_name,
    Total_Hubs = n_hubs,
    Differential_Hubs = n_diff_in_hubs,
    Risk_Hubs = n_risk_in_hubs,
    Protective_Hubs = n_protective_in_hubs,
    Percent_Differential = round(n_diff_in_hubs / n_hubs * 100, 2),
    stringsAsFactors = FALSE
  )

  if (n_diff_in_hubs > 0) {
    write.csv(hub_differential,
              file = paste0(outdir, outflag, group_name, "_differential_hubs.csv"),
              row.names = FALSE)
  }

  list(
    summary = result_df,
    details = hub_differential,
    risk_hubs = hub_differential$species_taxid[hub_differential$sig_marker == "risk"],
    protective_hubs = hub_differential$species_taxid[hub_differential$sig_marker == "protective"]
  )
}
