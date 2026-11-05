suppressMessages({library(ieugwasr); library(dplyr)})
args <- commandArgs(trailingOnly=TRUE)
start <- as.integer(args[1]); end <- as.integer(args[2]); budget <- as.numeric(args[3])
outdir <- file.path(Sys.getenv("NPM_MR_ROOT", "."), "data", "ivs")
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
ids <- paste0("ebi-a-GCST", start:end)
t0 <- Sys.time()
done <- 0
for (id in ids) {
  f <- file.path(outdir, paste0(id, ".csv"))
  if (file.exists(f)) next
  if (as.numeric(difftime(Sys.time(), t0, units="secs")) > budget) break
  ok <- FALSE
  for (attempt in 1:3) {
    r <- tryCatch(tophits(id, pval=5e-8), error=function(e) NULL)
    if (!is.null(r)) { ok <- TRUE; break }
    Sys.sleep(2)
  }
  if (!ok) { cat("FAIL", id, "\n", file=stderr()); next }
  if (nrow(r) == 0) {
    write.csv(data.frame(rsid=character(), p=numeric(), beta=numeric(), se=numeric(),
                         ea=character(), nea=character(), eaf=numeric(), n=numeric(), id=character()),
              f, row.names=FALSE)
  } else {
    r$id <- id
    write.csv(r, f, row.names=FALSE)
  }
  done <- done + 1
  Sys.sleep(0.2)
}
cat("processed", done, "files this run; total", length(list.files(outdir, pattern="csv$")), "\n")
