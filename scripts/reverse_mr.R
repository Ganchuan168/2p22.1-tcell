suppressMessages({library(ieugwasr); library(dplyr); library(data.table)})
ROOT <- Sys.getenv("NPM_MR_ROOT", ".")
args <- commandArgs(trailingOnly = TRUE)
start <- as.integer(args[1]); end <- as.integer(args[2]); budget <- as.numeric(args[3])

snps <- tophits("finn-b-N14_INFLAMMBREAST", pval = 5e-6)
snps <- snps[!duplicated(snps$rsid), ]
cl <- tryCatch(ld_clump(data.frame(rsid = snps$rsid, p = snps$p, id = "NPM"),
                        plink_bin = Sys.which("plink1.9"),
                        bfile = Sys.getenv("PLINK_1KG_EUR", file.path(ROOT, "data/ld/EUR")), r2 = 0.001, kb = 10000),
               error = function(e) snps)
snps <- snps[snps$rsid %in% cl$rsid, ]
cat("reverse instruments:", nrow(snps), "\n")
write.csv(snps, file.path(ROOT, "data/reverse_instruments.csv"), row.names = FALSE)

ids <- paste0("ebi-a-GCST", start:end)
outdir <- file.path(ROOT, "data/reverse")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
t0 <- Sys.time(); n <- 0
for (id in ids) {
  f <- file.path(outdir, paste0(id, ".csv"))
  if (file.exists(f)) next
  if (as.numeric(difftime(Sys.time(), t0, units = "secs")) > budget) break
  r <- tryCatch(associations(variants = snps$rsid, id = id), error = function(e) NULL)
  if (is.null(r)) next
  write.csv(r, f, row.names = FALSE)
  n <- n + 1
}
cat("done", n, "total", length(list.files(outdir, pattern = "csv$")), "\n")
