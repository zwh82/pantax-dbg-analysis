.libPaths(c("/home/wenhai/miniconda3/envs/r_4.3/lib/R/library","/home/wenhai/application/Qmatey/tools/R",.libPaths()))
suppressPackageStartupMessages({library(igraph); library(ggplot2)})
args <- commandArgs(trailingOnly=TRUE)
old_dir <- args[1]; new_dir <- args[2]
files_old <- sort(list.files(old_dir, recursive=TRUE))
files_new <- sort(list.files(new_dir, recursive=TRUE))
fail <- 0; checked <- 0
report <- function(ok, msg) { checked <<- checked+1; if(!ok){fail<<-fail+1; cat("FAIL:",msg,"\n")} }
if (!identical(files_old, files_new)) {
  cat("File set mismatch\n"); print(setdiff(files_old, files_new)); print(setdiff(files_new, files_old)); fail <- fail+1
}
for (f in intersect(files_old, files_new)) {
  if (grepl("[.]png$", f)) next
  a <- file.path(old_dir,f); b <- file.path(new_dir,f)
  if (grepl("[.]rds$", f)) {
    oa <- readRDS(a); ob <- readRDS(b)
    if (grepl("g_.*_full_network[.]rds$", f)) {
      ea <- as_data_frame(oa,what="edges"); eb <- as_data_frame(ob,what="edges")
      ea <- ea[order(ea$from,ea$to),]; eb <- eb[order(eb$from,eb$to),]
      va <- as_data_frame(oa,what="vertices"); vb <- as_data_frame(ob,what="vertices")
      va <- va[order(va$name),]; vb <- vb[order(vb$name),]
      ok <- identical(va,vb) && identical(ea,eb) && identical(graph_attr(oa),graph_attr(ob))
      report(ok, paste("graph",f))
    } else if (inherits(oa,"gg") && inherits(ob,"gg")) {
      ba <- tryCatch(ggplot_build(oa),error=function(e)NULL)
      bb <- tryCatch(ggplot_build(ob),error=function(e)NULL)
      ok <- isTRUE(all.equal(oa$labels,ob$labels)) && isTRUE(all.equal(oa$data,ob$data)) && isTRUE(all.equal(oa$mapping,ob$mapping))
      if (!is.null(ba) && !is.null(bb)) ok <- ok && isTRUE(all.equal(ba$data,bb$data))
      report(ok, paste("ggplot",f))
    } else {
      report(isTRUE(all.equal(oa,ob)), paste("rds",f))
    }
  } else if (grepl("[.](csv|tsv|txt)$", f)) {
    ta <- readLines(a, warn=FALSE); tb <- readLines(b, warn=FALSE)
    if (identical(ta,tb)) { report(TRUE, paste("text",f)); next }
    # numeric/table comparison with tolerance
    sep <- if (grepl("[.]tsv$",f)) "\t" else ","
    da <- tryCatch(read.table(a, sep=sep, header=TRUE, check.names=FALSE, stringsAsFactors=FALSE, comment.char="", quote="\""), error=function(e)NULL)
    db <- tryCatch(read.table(b, sep=sep, header=TRUE, check.names=FALSE, stringsAsFactors=FALSE, comment.char="", quote="\""), error=function(e)NULL)
    if (!is.null(da) && !is.null(db)) {
      report(isTRUE(all.equal(da,db,tolerance=1e-8,check.attributes=FALSE)), paste("table",f))
    } else {
      report(FALSE, paste("text",f))
    }
  }
}
cat(sprintf("\nChecked %d files, %d failures\n", checked, fail))
quit(status=if(fail>0) 1 else 0)
