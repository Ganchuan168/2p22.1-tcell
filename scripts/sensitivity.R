suppressMessages({library(data.table); library(dplyr); library(TwoSampleMR); library(MRPRESSO)})
ROOT <- Sys.getenv("NPM_MR_ROOT", ".")
source(file.path(ROOT, "scripts/mr_core.R"))
d <- load_data("outcome_N14.csv")
res <- read.csv(file.path(ROOT, "results/forward_N14.csv"))
hits <- res %>% filter(!is.na(p)) %>% arrange(p) %>% filter(p < 0.05)

out <- list()
for (i in seq_len(nrow(hits))) {
  t <- hits$id[i]
  m <- build_trait(d, t)
  if (is.null(m)) next
  r <- list(id = t, trait = hits$trait[i], nsnp = nrow(m))
  if (nrow(m) >= 3) {
    pr <- tryCatch({
      x <- mr_presso(BetaOutcome = "beta_out", BetaExposure = "beta.x", SdOutcome = "se_out",
                     SdExposure = "se", data = as.data.frame(m), OUTLIERtest = TRUE,
                     DISTORTIONtest = TRUE, NbDistribution = 1000, SignifThreshold = 0.05)
      list(global = x$`MR-PRESSO results`$`Global Test`$Pvalue,
           nout = length(x$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`))
    }, error = function(e) list(global = NA, nout = NA))
    r$presso_global_p <- pr$global; r$presso_outliers <- pr$nout
    loo <- tryCatch({
      dat <- data.frame(SNP = m$rsid, beta.exposure = m$beta.x, se.exposure = m$se,
                        beta.outcome = m$beta_out, se.outcome = m$se_out,
                        effect_allele.exposure = m$ea, other_allele.exposure = m$nea,
                        eaf.exposure = m$eaf, pval.exposure = m$p, id.exposure = m$id,
                        id.outcome = "out", outcome = "out", exposure = "exp", mr_keep = TRUE)
      l <- mr_leaveoneout(dat); max(abs(l$b - mean(l$b))) / sd(l$b)
    }, error = function(e) NA)
    r$loo_max_z <- loo
  } else { r$presso_global_p <- NA; r$presso_outliers <- NA; r$loo_max_z <- NA }
  # Steiger
  r2exp <- sum(m$R2, na.rm = TRUE)
  r2out <- sum(2 * m$af_alt * (1 - m$af_alt) * m$beta_out^2, na.rm = TRUE)
  st <- tryCatch(mr_steiger(m$p, m$pval, 3658, 235752, r2exp, r2out), error = function(e) NULL)
  r$steiger_p <- if (is.null(st)) NA else st$steiger_pval
  r$steiger_dir <- if (is.null(st)) NA else st$correct_causal_direction
  out[[i]] <- as.data.frame(r, stringsAsFactors = FALSE)
}
out <- bind_rows(out)
write.csv(out, file.path(ROOT, "results/sensitivity_hits.csv"), row.names = FALSE)
print(out, row.names = FALSE)
