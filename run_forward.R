source(file.path(Sys.getenv("NPM_MR_ROOT", "."), "scripts", "mr_core.R"))
args <- commandArgs(trailingOnly = TRUE)
outcome_file <- args[1]; tag <- args[2]
d <- load_data(outcome_file)
traits <- unique(d$iv$id)
one <- function(t) {
  m <- build_trait(d, t)
  r <- run_mr(m)
  if (is.null(r)) {
    data.frame(id = t, trait = unique(d$iv$trait[d$iv$id == t])[1],
               nsnp = 0, Fmean = NA, Fmin = NA, R2 = NA, method = NA,
               b = NA, se = NA, p = NA, stringsAsFactors = FALSE)
  } else {
    r$id <- t
    r$trait <- unique(d$iv$trait[d$iv$id == t])[1]
    r
  }
}
# mclapply is fork-based and unavailable on Windows; fall back to serial there.
run_traits <- function(x, fun, cores = 4) {
  if (.Platform$OS.type == "windows") {
    message("parallel::mclapply is unavailable on Windows; running serially")
    lapply(x, fun)
  } else {
    parallel::mclapply(x, fun, mc.cores = cores)
  }
}
res <- bind_rows(run_traits(traits, one))
res$OR <- exp(res$b); res$L95 <- exp(res$b - 1.96 * res$se); res$U95 <- exp(res$b + 1.96 * res$se)
res$FDR <- p.adjust(res$p, method = "BH")
res$Bonf <- res$p < 0.05 / sum(!is.na(res$p))
res$direction <- ifelse(res$b > 0, "risk", "protective")
write.csv(res, file.path(Sys.getenv("NPM_MR_ROOT", "."), "results", paste0("forward_", tag, ".csv")), row.names = FALSE)
cat("saved", nrow(res), "rows\n")
