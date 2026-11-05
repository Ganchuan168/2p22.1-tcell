suppressMessages({library(data.table); library(dplyr)})
ROOT <- Sys.getenv("NPM_MR_ROOT", ".")
comp <- c(A="T", T="A", C="G", G="C")
harm <- function(de, oe, do, oo, bo) {
  de <- toupper(de); oe <- toupper(oe); do <- toupper(do); oo <- toupper(oo)
  if (any(is.na(c(de,oe,do,oo))) || any(nchar(c(de,oe,do,oo)) != 1)) return(NA)
  if (de == do && oe == oo) return(bo)
  if (de == oo && oe == do) return(-bo)
  if (!is.na(comp[do]) && de == comp[do] && oe == comp[oo]) return(bo)
  if (!is.na(comp[oo]) && de == comp[oo] && oe == comp[do]) return(-bo)
  NA
}
iv <- fread(file.path(ROOT, "data/reverse_instruments.csv"))
iv$F <- (iv$beta / iv$se)^2
res <- list()
fs <- list.files(file.path(ROOT, "data/reverse"), full.names = TRUE)
for (f in fs) {
  o <- fread(f)
  # only directly measured associations are used; proxy rows returned by the
  # OpenGWAS query are discarded (the forward analysis likewise uses no proxies)
  n_proxy <- sum(o$proxy, na.rm = TRUE)
  o <- o[o$proxy == FALSE, ]
  if (nrow(o) == 0) next
  m <- merge(iv, o, by = "rsid")
  if (nrow(m) == 0) next
  m$b_out <- mapply(harm, m$ea.x, m$nea.x, m$ea.y, m$nea.y, m$beta.y)
  m <- m[!is.na(m$b_out) & !is.na(m$se.y) & m$se.y > 0, ]
  if (nrow(m) == 0) next
  b <- m$b_out / m$beta.x; s <- m$se.y / abs(m$beta.x)
  if (nrow(m) == 1) {
    est <- b; se <- s; p <- 2 * pnorm(-abs(est / se))
  } else {
    w <- 1 / s^2
    est <- sum(b * w) / sum(w); se <- sqrt(1 / sum(w)); p <- 2 * pnorm(-abs(est / se))
  }
  res[[length(res) + 1]] <- data.frame(id = o$id[1], trait = o$trait[1], nsnp = nrow(m), n_proxy_dropped = n_proxy,
                                       b = est, se = se, p = p, Fmean = mean(m$F),
                                       stringsAsFactors = FALSE)
}
res <- bind_rows(res)
res$OR <- exp(res$b); res$L95 <- exp(res$b - 1.96 * res$se); res$U95 <- exp(res$b + 1.96 * res$se)
res$FDR <- p.adjust(res$p, method = "BH")
res <- res[order(res$p), ]
write.csv(res, file.path(ROOT, "results/reverse_NPM.csv"), row.names = FALSE)
cat("traits:", nrow(res), "| p<0.05:", sum(res$p < 0.05), "| FDR<0.05:", sum(res$FDR < 0.05), "\n")
print(head(res[, c("trait", "nsnp", "OR", "L95", "U95", "p", "FDR")], 10), row.names = FALSE, digits = 3)
