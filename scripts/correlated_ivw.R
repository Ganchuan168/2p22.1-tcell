# Correlation-aware (GLS) IVW sensitivity analysis
suppressMessages({library(data.table); library(Matrix)})
ROOT <- Sys.getenv("NPM_MR_ROOT", ".")
s <- fread(file.path(ROOT, "tables/Supplementary_Table_S1_instrument_level.csv"))
bim <- fread(file.path(ROOT, "data/ld/all759.bim"), header = FALSE)
setnames(bim, c("chr", "rsid", "cm", "pos", "a1", "a2"))
R <- as.matrix(fread(file.path(ROOT, "data/ld/all759ld.ld"), header = FALSE))
rownames(R) <- bim$rsid; colnames(R) <- bim$rsid
f <- read.csv(file.path(ROOT, "results/forward_N14.csv")); f <- f[!is.na(f$p), ]

res <- list()
for (i in seq_len(nrow(f))) {
  g <- s[trait_id == f$id[i]]
  g <- g[rsid %in% rownames(R)]
  g <- g[!duplicated(g$rsid)]
  if (nrow(g) < 2) next
  Rg <- R[g$rsid, g$rsid, drop = FALSE]
  V <- diag(g$se_outcome) %*% Rg %*% diag(g$se_outcome)
  X <- matrix(g$beta_exposure, ncol = 1); Y <- matrix(g$beta_outcome, ncol = 1)
  Vi <- tryCatch(solve(V), error = function(e) NULL)
  if (is.null(Vi)) next
  b <- solve(t(X) %*% Vi %*% X) %*% (t(X) %*% Vi %*% Y)
  se <- sqrt(solve(t(X) %*% Vi %*% X)[1, 1])
  res[[length(res) + 1]] <- data.table(id = f$id[i], trait = f$trait[i],
    p_standard = f$p[i], b_standard = f$b[i], se_standard = f$se[i],
    b_gls = as.numeric(b), se_gls = as.numeric(se), k = nrow(g),
    max_r2 = max(Rg[upper.tri(Rg)]))
}
r <- rbindlist(res)
r[, p_gls := 2 * pnorm(-abs(b_gls / se_gls))]
r[, relative_change := abs(b_gls - b_standard) / abs(b_standard)]
fwrite(r, file.path(ROOT, "results/correlated_ivw_sensitivity.csv"))
cat("phenotypes:", nrow(r), " max pairwise r2:", max(r$max_r2), "\n")
cat("max relative change among FDR<0.2:",
    max(r[id %in% f$id[!is.na(f$FDR) & f$FDR < 0.2], relative_change]), "\n")
