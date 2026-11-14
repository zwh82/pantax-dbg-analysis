# Control-disease edge comparison and paired loss tests.

# Preserved requires the edge to still exist AND both endpoint abundance
# classes to be unchanged.  `edge_only_status` ignores the abundance check and
# is reported alongside as a sensitivity measure.
compare_pair <- function(control, disease) {
  dmatch <- match(control$edge_id, disease$edge_id)
  exists <- !is.na(dmatch)
  same_abundance <- rep(FALSE, nrow(control))
  same_abundance[exists] <-
    control$source_abundance[exists] == disease$source_abundance[dmatch[exists]] &
    control$target_abundance[exists] == disease$target_abundance[dmatch[exists]]
  control$status <- ifelse(exists & same_abundance, "preserved", "lost")
  control$edge_only_status <- ifelse(exists, "preserved", "lost")
  gained <- disease[!disease$edge_id %in% control$edge_id, , drop = FALSE]
  gained$status <- "gained"
  gained$edge_only_status <- "gained"
  list(control = control, gained = gained)
}

loss_summary <- function(control_status) {
  masks <- edge_masks(control_status)
  bind_rows_safe(purrr::imap(masks, function(mask, category) {
    x <- control_status[mask, , drop = FALSE]
    dplyr::tibble(
      edge_category = category,
      n_control_edges = nrow(x),
      n_preserved = sum(x$status == "preserved"),
      n_lost = sum(x$status == "lost"),
      loss_rate = safe_fraction(n_lost, n_control_edges),
      preservation_rate = safe_fraction(n_preserved, n_control_edges),
      edge_only_preservation_rate =
        safe_fraction(sum(x$edge_only_status == "preserved"), n_control_edges)
    )
  }))
}

status_core_fraction <- function(x) {
  x %>%
    dplyr::filter(edge_abundance_type == "abundant-rare") %>%
    dplyr::group_by(status) %>%
    dplyr::summarise(
      n_abundant_rare_edges = dplyr::n(),
      n_core_rare_abundant_edges = sum(is_core_rare_associated),
      fraction_core_rare = safe_fraction(
        n_core_rare_abundant_edges, n_abundant_rare_edges
      ),
      .groups = "drop"
    )
}

paired_loss_tests <- function(loss) {
  comparisons <- dplyr::tribble(
    ~group1, ~group2,
    "core_rare_associated", "noncore_rare_associated",
    "core_associated", "noncore_noncore",
    "core_rare_abundant", "other_abundant_rare",
    "core_abundant_associated", "noncore_abundant_associated"
  )
  out <- loss %>%
    dplyr::filter(edge_category %in% unique(c(comparisons$group1, comparisons$group2))) %>%
    dplyr::select(disease_name, cohort, edge_category, loss_rate) %>%
    tidyr::pivot_wider(names_from = edge_category, values_from = loss_rate) %>%
    dplyr::group_by(disease_name) %>%
    dplyr::group_modify(~ bind_rows_safe(purrr::map(seq_len(nrow(comparisons)), function(i) {
      a <- comparisons$group1[i]
      b <- comparisons$group2[i]
      keep <- complete.cases(.x[[a]], .x[[b]])
      if (sum(keep) < 2) {
        return(dplyr::tibble(group1 = a, group2 = b, n_pairs = sum(keep), p_value = NA_real_))
      }
      p <- wilcox.test(.x[[a]][keep], .x[[b]][keep],
        paired = TRUE, exact = FALSE
      )$p.value
      dplyr::tibble(group1 = a, group2 = b, n_pairs = sum(keep), p_value = p)
    }))) %>%
    dplyr::ungroup()
  out$p_adjust_bh <- p.adjust(out$p_value, method = "BH")
  out
}
