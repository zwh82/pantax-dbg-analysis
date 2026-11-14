# Per-group network figure and node table.

graph_plot <- function(g, outdir, group_name, abundance_colors, abundance_labels,
                       outflag = "") {
  node_table <- data.frame(
    Species = V(g)$name,
    AbundanceCategory = V(g)$abundance_category,
    Degree = V(g)$degree,
    Positive_Edges = V(g)$pos_edges,
    Negative_Edges = V(g)$neg_edges,
    Pos_Ratio = round(V(g)$pos_ratio, 3),
    IsHub = V(g)$is_hub,
    IsDifferential = V(g)$is_differential,
    DiffColor = ifelse(V(g)$is_differential,
                       ifelse(V(g)$label_color == "red", "Risk", "Protective"),
                       "None"),
    stringsAsFactors = FALSE
  )
  write.csv(node_table,
            file = paste0(outdir, outflag, group_name, "_node_info.csv"),
            row.names = FALSE)

  set.seed(123)
  layout <- layout_with_fr(g, weights = E(g)$abs_weight, niter = 2000)
  safe_outdir <- gsub("%", "%%", outdir)

  png_filename <- paste0(safe_outdir, outflag, group_name, "_network_abundance_colored.png")
  png(png_filename, width = 2800, height = 2400, res = 200)

  plot(g,
       layout = layout,
       vertex.color = V(g)$abundance_color,
       vertex.size = V(g)$size,
       vertex.frame.color = V(g)$hub_frame_color,
       vertex.frame.width = V(g)$hub_frame_width,
       vertex.label = ifelse(is.na(V(g)$label), "", V(g)$label),
       vertex.label.color = V(g)$label_color,
       vertex.label.cex = 0.6,
       vertex.label.dist = 0.5,
       vertex.label.font = 2,
       edge.width = E(g)$width,
       edge.color = E(g)$color,
       edge.curved = 0.2,
       main = paste0(group_name, " Network by abundance\n(Stability: ",
                     round(g$Stability, 4), ")"))

  legend("topright",
         legend = c(abundance_labels[["high"]], abundance_labels[["low"]],
                    abundance_labels[["very_low"]], abundance_labels[["other"]],
                    "Hub nodes"),
         col = c(abundance_colors["high"], abundance_colors["low"],
                 abundance_colors["very_low"], abundance_colors["other"], "gold"),
         pch = 19,
         pt.cex = 1.5,
         bty = "n",
         title = "Abundance Category & Hubs")

  legend("bottomright",
         legend = c("Positive correlation", "Negative correlation"),
         col = c("#f783ac", "#3bc9db"),
         lwd = 3,
         bty = "n")

  dev.off()

  saveRDS(g, paste0(outdir, "g_", group_name, "_full_network.rds"))

  cat(group_name, "network processing completed.\n")
  cat("Nodes:", vcount(g), "Edges:", ecount(g), "\n")
  cat("Hub nodes:", sum(V(g)$is_hub), "\n")
  cat("Differential hubs:", sum(V(g)$is_hub & V(g)$is_differential), "\n")
  cat("  Risk hubs:", sum(V(g)$is_hub & V(g)$is_differential & V(g)$label_color == "red"), "\n")
  cat("  Protective hubs:", sum(V(g)$is_hub & V(g)$is_differential & V(g)$label_color == "green"), "\n")
  cat("Output files saved in:", outdir, "\n")
}
