# Variant-level statistics for the six core variants at 2p22.1
suppressMessages(library(data.table))
ROOT <- Sys.getenv("NPM_MR_ROOT", ".")
s <- fread(file.path(ROOT, "tables/Supplementary_Table_S1_instrument_level.csv"))
s[, wald := beta_outcome / beta_exposure]
s[, wald_se := se_outcome / abs(beta_exposure)]
s[, wp := 2 * pnorm(-abs(wald / wald_se))]
v <- s[, .(minp = min(wp)), by = rsid][, fdr := p.adjust(minp, "BH")]

bim <- fread(file.path(ROOT, "data/ld/all759.bim"), header = FALSE)
setnames(bim, c("chr", "rsid", "cm", "pos", "a1", "a2"))
R <- as.matrix(fread(file.path(ROOT, "data/ld/all759ld.ld"), header = FALSE))
rownames(R) <- bim$rsid; colnames(R) <- bim$rsid

core <- c("rs6707421", "rs12712610", "rs6751481", "rs413431", "rs408686", "rs17607399")
out <- rbindlist(lapply(core, function(id) {
  x <- s[rsid == id][which.min(wp)]
  data.table(variant = id, cluster = c("A1","A1","A2","B","B","C")[match(id, core)],
    position_grch38 = x$position_grch38, ea = x$effect_allele, oa = x$other_allele,
    eaf = x$eaf, b_exp = x$beta_exposure, se_exp = x$se_exposure, p_exp = x$p_exposure,
    b_out = x$beta_outcome, se_out = x$se_outcome, OR = exp(x$wald),
    L95 = exp(x$wald - 1.96 * x$wald_se), U95 = exp(x$wald + 1.96 * x$wald_se),
    p = x$wp, fdr = v[rsid == id, fdr], phenotype = x$phenotype,
    r2_lead = ifelse(id == "rs6707421", 1, R[id, "rs6707421"]))
}))
fwrite(out, file.path(ROOT, "results/variant_level_table.csv"))
print(out, digits = 3)
