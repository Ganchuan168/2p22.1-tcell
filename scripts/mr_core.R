suppressMessages({
  library(data.table); library(dplyr); library(TwoSampleMR)
})

ROOT <- Sys.getenv("NPM_MR_ROOT", ".")
comp <- c(A="T", T="A", C="G", G="C")

harmonise_one <- function(de, oe, do, oo, bo, so, fo) {
  # de/oe exposure effect/other allele; do/oo outcome; bo/so beta/se for outcome alt; fo af of outcome alt
  if (is.na(de) || is.na(oe) || is.na(do) || is.na(oo)) return(NULL)
  de <- toupper(de); oe <- toupper(oe); do <- toupper(do); oo <- toupper(oo)
  if (nchar(de) != 1 || nchar(oe) != 1 || nchar(do) != 1 || nchar(oo) != 1) return(NULL)
  pal <- (de %in% c("A","T") && oe %in% c("A","T")) || (de %in% c("C","G") && oe %in% c("C","G"))
  if (de == do && oe == oo) {
    if (pal && !is.na(fo) && fo > 0.42 && fo < 0.58) return(NULL)
    return(c(beta = bo, se = so))
  }
  if (de == oo && oe == do) {
    if (pal && !is.na(fo) && fo > 0.42 && fo < 0.58) return(NULL)
    return(c(beta = -bo, se = so))
  }
  # try complement of outcome alleles
  do2 <- comp[do]; oo2 <- comp[oo]
  if (is.na(do2) || is.na(oo2)) return(NULL)
  if (de == do2 && oe == oo2) {
    fo2 <- ifelse(is.na(fo), NA, 1 - fo)
    if (pal && !is.na(fo2) && fo2 > 0.42 && fo2 < 0.58) return(NULL)
    return(c(beta = bo, se = so))
  }
  if (de == oo2 && oe == do2) {
    fo2 <- ifelse(is.na(fo), NA, 1 - fo)
    if (pal && !is.na(fo2) && fo2 > 0.42 && fo2 < 0.58) return(NULL)
    return(c(beta = -bo, se = so))
  }
  NULL
}

load_data <- function(outcome_file) {
  iv <- rbindlist(lapply(list.files(file.path(ROOT, "data/ivs"), full.names = TRUE), fread),
                  fill = TRUE) %>% as.data.frame()
  out <- fread(file.path(ROOT, "data", outcome_file)) %>% as.data.frame()
  out <- out[!duplicated(out$query_rsid), ]
  list(iv = iv, out = out)
}

build_trait <- function(d, trait_id) {
  x <- d$iv[d$iv$id == trait_id, ]
  if (nrow(x) == 0) return(NULL)
  x <- x[!is.na(x$eaf) & x$eaf >= 0.01 & x$eaf <= 0.99, ]
  if (nrow(x) == 0) return(NULL)
  if (nrow(x) > 1) {
    cl <- tryCatch(ld_clump(data.frame(rsid = x$rsid, p = x$p, id = x$id),
                            plink_bin = Sys.which("plink1.9"),
                            bfile = Sys.getenv("PLINK_1KG_EUR", file.path(ROOT, "data/ld/EUR")),
                            r2 = 0.001, kb = 10000), error = function(e) NULL)
    if (!is.null(cl)) x <- x[x$rsid %in% cl$rsid, ]
  }
  if (nrow(x) == 0) return(NULL)
  m <- merge(x, d$out, by.x = "rsid", by.y = "query_rsid", all.x = FALSE)
  if (nrow(m) == 0) return(NULL)
  res <- t(sapply(seq_len(nrow(m)), function(i) {
    h <- harmonise_one(m$ea[i], m$nea[i], m$alt[i], m$ref[i], m$beta.y[i], m$sebeta[i], m$af_alt[i])
    if (is.null(h)) return(c(beta = NA_real_, se = NA_real_))
    as.numeric(h[c("beta", "se")])
  }))
  m$beta_out <- res[, 1]; m$se_out <- res[, 2]
  m <- m[!is.na(m$beta_out) & !is.na(m$se_out) & m$se_out > 0, ]
  if (nrow(m) == 0) return(NULL)
  m$F <- (m$beta.x / m$se)^2
  m$R2 <- 2 * m$eaf * (1 - m$eaf) * m$beta.x^2
  m
}

run_mr <- function(m) {
  if (is.null(m) || nrow(m) < 1) return(NULL)
  nsnp <- nrow(m)
  out <- list(nsnp = nsnp, Fmean = mean(m$F, na.rm = TRUE), Fmin = min(m$F, na.rm = TRUE),
              R2 = sum(m$R2, na.rm = TRUE))
  if (nsnp == 1) {
    b <- m$beta_out / m$beta.x
    s <- m$se_out / abs(m$beta.x)
    out$method <- "Wald ratio"; out$b <- b; out$se <- s; out$p <- 2 * pnorm(-abs(b / s))
  } else {
    m <- m[complete.cases(m[,c("beta.x","se","beta_out","se_out","eaf")]), ]
    if (nrow(m) < 2) return(NULL)
    dat <- data.frame(SNP = m$rsid, beta.exposure = m$beta.x, se.exposure = m$se,
                      beta.outcome = m$beta_out, se.outcome = m$se_out,
                      effect_allele.exposure = m$ea, other_allele.exposure = m$nea,
                      eaf.exposure = m$eaf, pval.exposure = m$p,
                      id.exposure = m$id, id.outcome = "out", outcome = "out", exposure = "exp",
                      mr_keep = TRUE)
    b_exp <- m$beta.x; b_out <- m$beta_out; se_exp <- m$se; se_out <- m$se_out
    ivw <- tryCatch(mr_ivw(b_exp, b_out, se_exp, se_out), error = function(e) NULL)
    egg <- tryCatch(mr_egger_regression(b_exp, b_out, se_exp, se_out), error = function(e) NULL)
    wm <- tryCatch(mr_weighted_median(b_exp, b_out, se_exp, se_out), error = function(e) NULL)
    wmode <- tryCatch(mr_weighted_mode(b_exp, b_out, se_exp, se_out), error = function(e) NULL)
    smode <- tryCatch(mr_simple_mode(b_exp, b_out, se_exp, se_out), error = function(e) NULL)
    if (is.null(ivw)) return(NULL)
    out$method <- "IVW"; out$b <- ivw$b; out$se <- ivw$se; out$p <- ivw$pval
    out$egger_b <- if (is.null(egg)) NA else egg$b
    out$egger_p <- if (is.null(egg)) NA else egg$pval
    out$wm_b <- if (is.null(wm)) NA else as.numeric(wm$b)
    out$wm_p <- if (is.null(wm)) NA else as.numeric(wm$pval)
    out$wmode_b <- if (is.null(wmode)) NA else as.numeric(wmode$b)
    out$wmode_p <- if (is.null(wmode)) NA else as.numeric(wmode$pval)
    out$smode_b <- if (is.null(smode)) NA else as.numeric(smode$b)
    out$smode_p <- if (is.null(smode)) NA else as.numeric(smode$pval)
    q <- tryCatch(mr_heterogeneity(dat, method_list = "mr_ivw"), error = function(e) NULL)
    out$Q_p <- if (is.null(q) || nrow(q) == 0) NA else q$Q_pval[1]
    pl <- tryCatch(mr_pleiotropy_test(dat), error = function(e) NULL)
    out$egger_intercept <- if (is.null(pl)) NA else pl$egger_intercept
    out$egger_intercept_p <- if (is.null(pl)) NA else pl$pval
  }
  as.data.frame(out, stringsAsFactors = FALSE)
}
