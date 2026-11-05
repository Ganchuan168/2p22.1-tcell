suppressMessages({library(data.table); library(dplyr); library(openxlsx)})
ROOT <- Sys.getenv("NPM_MR_ROOT", ".")
source(file.path(ROOT, "scripts/mr_core.R"))
d <- load_data("outcome_N14.csv")
f <- read.csv(file.path(ROOT, "results/forward_N14.csv"))
o <- read.csv(file.path(ROOT, "results/forward_O15.csv"))
h <- f[!is.na(f$FDR) & f$FDR < 0.2, ]
h <- h[order(h$FDR, h$p), ]
lead <- sapply(h$id, function(t) { m <- build_trait(d, t); paste(m$rsid, collapse = "; ") })
h$lead_snp <- lead
h <- merge(h, o[, c("id", "p", "FDR", "b", "se")], by = "id", suffixes = c("", "_O15"))
h <- h[order(h$FDR, h$p), ]
# difference between NPM and lactational mastitis estimates (conservative: assumes independence)
h$z_diff <- (h$b - h$b_O15) / sqrt(h$se^2 + h$se_O15^2)
h$p_diff <- 2 * pnorm(-abs(h$z_diff))

tab2 <- data.frame(
  Phenotype = h$trait,
  `No. of IVs` = h$nsnp,
  `Lead variant` = h$lead_snp,
  OR = sprintf("%.2f", h$OR), `95% CI` = sprintf("%.2f-%.2f", h$L95, h$U95),
  P = signif(h$p, 3), FDR = signif(h$FDR, 3),
  Bonferroni = ifelse(h$Bonf, "Yes", "No"),
  `P, lactational mastitis` = signif(h$p_O15, 3),
  `P for difference` = signif(h$p_diff, 3),
  `Mean F` = round(h$Fmean, 1),
  check.names = FALSE)

rev <- read.csv(file.path(ROOT, "results/reverse_NPM.csv"))
tab3 <- head(rev, 15)[, c("trait", "nsnp", "OR", "L95", "U95", "p", "FDR", "Fmean")]
names(tab3) <- c("Phenotype", "No. of IVs", "OR", "L95", "U95", "P", "FDR", "Mean F")
tab3[, 3:8] <- lapply(tab3[, 3:8], function(x) signif(x, 3))

tab1 <- data.frame(
  Characteristic = c("Immune cell phenotypes (exposure)", "Source", "Sample size",
                     "Phenotype classes", "NPM (primary outcome)", "Source", "Cases / controls",
                     "Definition", "Lactational mastitis (control outcome)", "Cases / controls",
                     "Definition", "Ancestry", "IV threshold", "LD clumping"),
  Value = c("731", "Orru et al., Nat Genet 2020 (GCST90001391-GCST90002121)", "3,757 Sardinian individuals",
            "AC 118; MFI 389; RC 192; MP 32", "Non-puerperal mastitis", "FinnGen R11",
            "2,167 / 233,585", "ICD-10 N61 (inflammatory disorders of breast)",
            "Infections of breast associated with childbirth", "1,233 / 232,519",
            "ICD-10 O91", "European", "P < 5 x 10^-8", "r2 < 0.001, 10,000 kb"),
  check.names = FALSE)

wb <- createWorkbook()
addWorksheet(wb, "Table1"); writeData(wb, "Table1", tab1)
addWorksheet(wb, "Table2"); writeData(wb, "Table2", tab2)
addWorksheet(wb, "Table3"); writeData(wb, "Table3", tab3)
addWorksheet(wb, "TableS1_forward"); writeData(wb, "TableS1_forward", f)
addWorksheet(wb, "TableS2_reverse"); writeData(wb, "TableS2_reverse", rev)
st <- createStyle(textDecoration = "bold")
addStyle(wb, "Table1", st, rows = 1, cols = 1:2, gridExpand = TRUE)
addStyle(wb, "Table2", st, rows = 1, cols = 1:ncol(tab2), gridExpand = TRUE)
addStyle(wb, "Table3", st, rows = 1, cols = 1:ncol(tab3), gridExpand = TRUE)
saveWorkbook(wb, file.path(ROOT, "results/tables.xlsx"), overwrite = TRUE)
write.csv(tab2, file.path(ROOT, "results/table2.csv"), row.names = FALSE)
cat("tables written; FDR<0.2 traits:", nrow(tab2), "\n")
print(tab2[, c("Phenotype", "OR", "95% CI", "P", "FDR", "Bonferroni")], row.names = FALSE)
