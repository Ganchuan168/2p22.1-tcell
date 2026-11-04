suppressMessages({library(ggplot2); library(dplyr); library(data.table); library(scales); library(ggrepel)})
ROOT <- Sys.getenv("NPM_MR_ROOT", ".")
f <- read.csv(file.path(ROOT, "results/forward_N14.csv"))
allf <- f

disp <- c(
  "Terminally Differentiated CD8+ T cell Absolute Count" = "Terminally differentiated CD8+ T cells, absolute count",
  "Terminally Differentiated CD8+ T cell %T cell" = "Terminally differentiated CD8+ T cells, % of T cells",
  "CD45RA- CD4+ T cell %T cell" = "CD45RA− CD4+ T cells, % of T cells",
  "Terminally Differentiated CD4-CD8- T cell %CD4-CD8- T cell" = "Terminally differentiated CD4−CD8− T cells, % of DN T cells",
  "CD28+ CD45RA+ CD8dim T cell %CD8dim T cell" = "CD28+CD45RA+ CD8dim T cells, % of CD8dim T cells",
  "CD28+ CD45RA+ CD8dim T cell Absolute Count" = "CD28+CD45RA+ CD8dim T cells, absolute count",
  "Activated & resting CD4 regulatory T cell %CD4 regulatory T cell" = "Activated and resting CD4 Tregs, % of Tregs",
  "Secreting CD4 regulatory T cell %CD4 regulatory T cell" = "Secreting CD4 Tregs, % of Tregs",
  "Activated & resting CD4 regulatory T cell %CD4+ T cell" = "Activated and resting CD4 Tregs, % of CD4+ T cells",
  "Resting CD4 regulatory T cell %CD4+ T cell" = "Resting CD4 Tregs, % of CD4+ T cells",
  "Activated & resting CD4 regulatory T cell Absolute Count" = "Activated and resting CD4 Tregs, absolute count",
  "Central Memory CD4+ T cell %CD4+ T cell" = "Central memory CD4+ T cells, % of CD4+ T cells",
  "HLA DR on CD14- CD16+ monocyte" = "HLA-DR on CD14−CD16+ monocytes",
  "FSC-A on CD14+ monocyte" = "FSC-A on CD14+ monocytes",
  "FSC-A on CD4+ T cell" = "FSC-A on CD4+ T cells",
  "HLA DR+ CD4+ T cell %lymphocyte" = "HLA-DR+ CD4+ T cells, % of lymphocytes",
  "CD62L on granulocyte" = "CD62L on granulocytes",
  "FSC-A on T cell" = "FSC-A on T cells",
  "SSC-A on CD4+ T cell" = "SSC-A on CD4+ T cells",
  "CCR7 on naive CD8+ T cell" = "CCR7 on naive CD8+ T cells",
  "Naive CD8+ T cell Absolute Count" = "Naive CD8+ T cells, absolute count",
  "HLA DR+ CD4+ T cell %T cell" = "HLA-DR+ CD4+ T cells, % of T cells",
  "Activated CD4 regulatory T cell %CD4 regulatory T cell" = "Activated CD4 Tregs, % of Tregs",
  "HLA DR+ CD4+ T cell Absolute Count" = "HLA-DR+ CD4+ T cells, absolute count",
  "Naive CD8+ T cell %T cell" = "Naive CD8+ T cells, % of T cells",
  "CD28+ CD45RA+ CD8+ T cell Absolute Count" = "CD28+CD45RA+ CD8+ T cells, absolute count",
  "SSC-A on T cell" = "SSC-A on T cells",
  "Activated CD4 regulatory T cell %CD4+ T cell" = "Activated CD4 Tregs, % of CD4+ T cells",
  "CD28+ CD45RA+ CD8+ T cell %T cell" = "CD28+CD45RA+ CD8+ T cells, % of T cells"
)

theme_pub <- theme_classic(base_size = 10) +
  theme(axis.text = element_text(colour = "black"),
        plot.title = element_text(face = "bold", size = 11),
        legend.position = "top")

# Figure 2 volcano
f2 <- f[!is.na(f$p), ]
f2$sig <- ifelse(!is.na(f2$FDR) & f2$FDR < 0.05, "FDR < 0.05",
           ifelse(!is.na(f2$FDR) & f2$FDR < 0.2, "0.05 ≤ FDR < 0.2",
           ifelse(f2$p < 0.05, "P < 0.05", "Not significant")))
f2$sig <- factor(f2$sig, levels = c("FDR < 0.05", "0.05 ≤ FDR < 0.2", "P < 0.05", "Not significant"))
keep_lab <- c("Terminally Differentiated CD8+ T cell %T cell",
              "CD45RA- CD4+ T cell %T cell",
              "Activated & resting CD4 regulatory T cell %CD4 regulatory T cell",
              "Secreting CD4 regulatory T cell %CD4 regulatory T cell",
              "Resting CD4 regulatory T cell %CD4+ T cell",
              "Central Memory CD4+ T cell %CD4+ T cell")
lab <- f2[f2$trait %in% keep_lab, ]
short_map <- c(
  "Terminally Differentiated CD8+ T cell Absolute Count" = "TEMRA CD8+ (count)",
  "Terminally Differentiated CD8+ T cell %T cell" = "TEMRA CD8+ (%T)",
  "CD45RA- CD4+ T cell %T cell" = "CD45RA− CD4+ (%T)",
  "Terminally Differentiated CD4-CD8- T cell %CD4-CD8- T cell" = "CD4−CD8− T (%DN)",
  "CD28+ CD45RA+ CD8dim T cell %CD8dim T cell" = "CD28+CD45RA+ CD8dim (%)",
  "Activated & resting CD4 regulatory T cell %CD4 regulatory T cell" = "Treg, act+rest (%Treg)",
  "Secreting CD4 regulatory T cell %CD4 regulatory T cell" = "Treg, secreting (%Treg)",
  "Resting CD4 regulatory T cell %CD4+ T cell" = "Treg, resting (%CD4)",
  "Central Memory CD4+ T cell %CD4+ T cell" = "Central memory CD4+ (%CD4)"
)
lab$short <- ifelse(lab$trait %in% names(short_map), short_map[lab$trait], substr(lab$trait, 1, 26))

p2 <- ggplot(f2, aes(b, -log10(p), colour = sig)) +
  geom_point(size = 1.4, alpha = 0.8) +
  scale_colour_manual(values = c("FDR < 0.05" = "#C0392B", "0.05 ≤ FDR < 0.2" = "#E67E22",
                                 "P < 0.05" = "#F1C40F", "Not significant" = "grey75")) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, colour = "grey40") +
  geom_hline(yintercept = -log10(0.05 / 578), linetype = 3, colour = "grey40") +
  geom_text_repel(data = lab, aes(b, -log10(p), label = short), size = 2.5, colour = "grey15",
                  min.segment.length = 0, segment.size = 0.2, segment.colour = "grey50",
                  box.padding = 0.45, max.overlaps = Inf, seed = 1) +
  labs(x = "MR estimate (log OR per SD of immune trait)", y = expression(-log[10](P)), colour = NULL,
       title = "Immune cell phenotypes and non-puerperal mastitis") +
  theme_pub
ggsave(file.path(ROOT, "figures/fig2_volcano.pdf"), p2, width = 6.9, height = 4.9, device = cairo_pdf)
ggsave(file.path(ROOT, "figures/fig2_volcano.png"), p2, width = 6.9, height = 4.9, dpi = 400, device = ragg::agg_png)

# Figure 3 forest
g <- allf[!is.na(allf$FDR) & allf$FDR < 0.2, ]
g <- g[order(g$FDR, g$p), ]
g$tier <- ifelse(g$FDR < 0.05, "FDR < 0.05", "0.05 ≤ FDR < 0.2")
g$tier <- factor(g$tier, levels = c("FDR < 0.05", "0.05 ≤ FDR < 0.2"))
g$trait_disp <- unname(disp[g$trait])
g$trait_disp <- factor(g$trait_disp, levels = rev(g$trait_disp))
p3 <- ggplot(g, aes(x = OR, y = trait_disp, colour = tier)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
  geom_errorbarh(aes(xmin = L95, xmax = U95), height = 0.2) +
  geom_point(size = 2) +
  scale_x_log10() +
  scale_colour_manual(values = c("FDR < 0.05" = "#C0392B", "0.05 ≤ FDR < 0.2" = "#E67E22")) +
  labs(x = "Odds ratio (95% CI) per SD increase", y = NULL, colour = NULL,
       title = "Immune phenotypes associated with NPM") +
  theme_pub
ggsave(file.path(ROOT, "figures/fig3_forest.pdf"), p3, width = 6.9, height = 4.6, device = cairo_pdf)
ggsave(file.path(ROOT, "figures/fig3_forest.png"), p3, width = 6.9, height = 4.6, dpi = 400, device = ragg::agg_png)
cat("fig2-3 done\n")
