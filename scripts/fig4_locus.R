suppressMessages({library(data.table); library(dplyr); library(ggplot2); library(ggrepel)})
ROOT <- Sys.getenv("NPM_MR_ROOT", ".")
source(file.path(ROOT, "scripts/mr_core.R"))
d <- load_data("outcome_N14.csv")
pos <- fread(file.path(ROOT, "data/positions_grch38.csv"))
f <- read.csv(file.path(ROOT, "results/forward_N14.csv"))
h <- f[!is.na(f$FDR) & f$FDR < 0.2, ]
h <- h[order(h$FDR, h$p), ]

rows <- list()
for (i in seq_len(nrow(h))) {
  m <- build_trait(d, h$id[i])
  if (is.null(m)) next
  m$wald <- m$beta_out / m$beta.x
  m$wald_se <- m$se_out / abs(m$beta.x)
  m$OR <- exp(m$wald)
  m$trait <- h$trait[i]
  rows[[i]] <- m[, c("rsid", "trait", "OR", "wald_se")]
}
w <- bind_rows(rows) %>% distinct(rsid, trait, OR, wald_se) %>% left_join(pos, by = "rsid")
w <- w[!is.na(w$pos) & w$chr == 2, ]
w$mb <- w$pos / 1e6
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
  "Central Memory CD4+ T cell %CD4+ T cell" = "Central memory CD4+ T cells, % of CD4+ T cells"
)
w$trait_disp <- unname(disp[w$trait])
w$trait_disp <- factor(w$trait_disp, levels = rev(unique(unname(disp[h$trait]))))
w$sig <- ifelse(w$OR > 1, "OR > 1", "OR < 1")

lead <- w %>% mutate(cluster = round(mb, 2)) %>% group_by(cluster) %>%
  summarise(mb = mean(mb), label = paste(sort(unique(rsid)), collapse = " / "), .groups = "drop") %>%
  arrange(mb)
lead$tag <- letters[seq_len(nrow(lead))]

p <- ggplot(w, aes(x = mb, y = trait_disp)) +
  geom_vline(data = lead, aes(xintercept = mb), colour = "grey80", linewidth = 0.3, linetype = 3) +
  geom_point(aes(colour = sig, shape = sig, size = 1 / wald_se), alpha = 0.9) +
  scale_shape_manual(values = c("OR > 1" = 17, "OR < 1" = 16), name = NULL) +
  geom_text(data = lead, aes(x = mb, y = Inf, label = tag), inherit.aes = FALSE,
            vjust = -0.6, size = 3.4, fontface = "bold", colour = "grey25") +
  scale_colour_manual(values = c("OR > 1" = "#C0392B", "OR < 1" = "#2471A3"), name = NULL) +
  guides(shape = guide_legend(override.aes = list(size = 2.2)), colour = guide_legend(override.aes = list(size = 2.2))) +
  scale_size_continuous(range = c(0.6, 2.8), guide = "none") +
  scale_x_continuous(limits = c(38.5, 42.0), breaks = seq(38.5, 42, 0.5),
                     labels = function(x) paste0(x, " Mb")) +
  coord_cartesian(clip = "off") +
  labs(x = "Chromosome 2 position (GRCh38)", y = NULL,
       title = "Twelve of the 14 phenotypes with FDR < 0.2\nmap to a 3-Mb interval at 2p22.1\u2013p21") +
  theme_classic(base_size = 9) +
  theme(axis.text = element_text(colour = "black"),
        axis.text.y = element_text(size = 7.2),
        plot.title = element_text(face = "bold", size = 8.2, hjust = 0, lineheight = 1.05),
        legend.position = "top",
        plot.margin = margin(26, 6, 4, 4))

ggsave(file.path(ROOT, "figures/fig4_locus.pdf"), p, width = 6.9, height = 4.3, device = cairo_pdf)
ggsave(file.path(ROOT, "figures/fig4_locus.png"), p, width = 6.9, height = 4.3, dpi = 400, device = ragg::agg_png)
write.csv(w, file.path(ROOT, "results/per_variant_wald.csv"), row.names = FALSE)
cat("fig4 done;", nrow(lead), "clusters:", paste(lead$label, collapse = " | "), "\n")
