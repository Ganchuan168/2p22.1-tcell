suppressMessages({library(ggplot2)})
ROOT <- Sys.getenv("NPM_MR_ROOT", ".")

mk <- function(x, y, w, h, label, fill) data.frame(x, y, w, h, label, fill)
b <- rbind(
  mk(0.12, 0.84, 0.215, 0.24, "Exposure\n731 circulating immune\ncell phenotypes\n3,757 Sardinian individuals", "#EAF2F8"),
  mk(0.36, 0.84, 0.215, 0.24, "Instrument selection\nP < 5×10⁻⁸; r² < 0.001\nmean F > 10\n612 phenotypes ≥1 variant", "#EAF2F8"),
  mk(0.60, 0.84, 0.215, 0.24, "Harmonisation\n578 phenotypes analysed\n1–25 variants\nmedian 2", "#F4F6F6"),
  mk(0.84, 0.84, 0.215, 0.24, "Forward MR\nIVW, MR-Egger, weighted\nmedian and mode,\nsimple mode", "#E8F8F5"),
  mk(0.60, 0.47, 0.215, 0.18, "Primary outcome\nNon-puerperal mastitis\nICD-10 N61\n2,167 / 233,585", "#FDEDEC"),
  mk(0.84, 0.47, 0.215, 0.18, "Comparator outcome\nLactational mastitis\nICD-10 O91\n1,233 / 232,519", "#FDEDEC"),
  mk(0.24, 0.13, 0.44, 0.20, "Reverse-direction MR\n7 NPM-associated variants (P < 5×10⁻⁶)\n1-4 per phenotype; 730 phenotypes tested", "#E8F8F5")
)

h_arrow <- data.frame(x = c(0.2275, 0.4675, 0.7075), y = 0.84,
                      xend = c(0.2525, 0.4925, 0.7325), yend = 0.84)
v_arrow <- data.frame(x = c(0.60, 0.84), y = c(0.56, 0.56), xend = c(0.60, 0.84), yend = c(0.74, 0.74))
d_arrow <- data.frame(x = 0.50, y = 0.38, xend = 0.42, yend = 0.21)

p <- ggplot() +
  geom_rect(data = b, aes(xmin = x - w/2, xmax = x + w/2, ymin = y - h/2, ymax = y + h/2, fill = fill),
            colour = "grey40", linewidth = 0.3) +
  geom_text(data = b, aes(x = x, y = y, label = label), size = 2.55, lineheight = 0.90, colour = "grey10") +
  geom_segment(data = h_arrow, aes(x = x, y = y, xend = xend, yend = yend),
               arrow = arrow(length = unit(0.08, "cm"), type = "closed"), linewidth = 0.3, colour = "grey30") +
  geom_segment(data = v_arrow, aes(x = x, y = y, xend = xend, yend = yend),
               arrow = arrow(length = unit(0.08, "cm"), type = "closed"), linewidth = 0.3, colour = "grey30") +
  geom_segment(data = d_arrow, aes(x = x, y = y, xend = xend, yend = yend), linetype = 2,
               arrow = arrow(length = unit(0.08, "cm"), type = "closed"), linewidth = 0.3, colour = "grey30") +
  scale_fill_identity() +
  coord_cartesian(xlim = c(0, 1), ylim = c(0.02, 1), expand = FALSE) +
  theme_void() + theme(plot.margin = margin(3, 3, 3, 3))

ggsave(file.path(ROOT, "figures/fig1_design.pdf"), p, width = 6.9, height = 3.6, device = cairo_pdf)
ggsave(file.path(ROOT, "figures/fig1_design.png"), p, width = 6.9, height = 3.6, dpi = 400, device = ragg::agg_png)
cat("fig1 done\n")
