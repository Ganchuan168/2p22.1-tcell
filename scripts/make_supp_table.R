suppressMessages({library(data.table); library(dplyr)})
ROOT <- Sys.getenv("NPM_MR_ROOT", ".")
source(file.path(ROOT, "scripts/mr_core.R"))

d <- load_data("outcome_N14.csv")
pos <- fread(file.path(ROOT, "data/positions_grch38.csv"))
pos <- pos[!duplicated(rsid), .(rsid, chr_grch38 = chr, position_grch38 = pos)]
f <- read.csv(file.path(ROOT, "results/forward_N14.csv"))
ids <- f$id[!is.na(f$p)]

res <- vector("list", length(ids))
for (k in seq_along(ids)) {
  m <- build_trait(d, ids[k])
  if (is.null(m)) next
  m$trait_name <- f$trait[match(ids[k], f$id)]
  m <- merge(m, pos, by = "rsid", all.x = TRUE, sort = FALSE)
  res[[k]] <- data.frame(
    phenotype = m$trait_name,
    trait_id = m$id,
    rsid = m$rsid,
    chr = m$chr_grch38,
    position_grch38 = m$position_grch38,
    effect_allele = m$ea,
    other_allele = m$nea,
    eaf = m$eaf,
    beta_exposure = m$beta.x,
    se_exposure = m$se,
    p_exposure = m$p,
    beta_outcome = m$beta_out,
    se_outcome = m$se_out,
    F_statistic = m$F,
    R2 = m$R2,
    stringsAsFactors = FALSE)
}
tab <- bind_rows(res) %>% arrange(trait_id, chr, position_grch38)
fwrite(tab, file.path(ROOT, "tables/Supplementary_Table_S1_instrument_level.csv"))
cat("rows", nrow(tab), "traits", length(unique(tab$trait_id)),
    "rsids", length(unique(tab$rsid)), "\n")
