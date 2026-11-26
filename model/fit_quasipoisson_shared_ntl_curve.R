## 23_klog60_driver.R — full-scope screen of the KLOG specification (60 configurations x 2 sexes).
## One process per (scope, architecture, period): raw cache read once per sex, k = 4..9 fitted in-process.
## Specification = production V2 expanded, except L enters as s(log L, k, bs="cr") with production knot rule
## applied on the log scale (city_deviation sz term likewise). Response, encoded zeros, rows, weights (rep 0),
## controls, family quasipoisson(log), fREML, discrete=10000, epsilon 1e-7, maxit 200: unchanged.
## Usage: Rscript 23_klog60_driver.R <scope> <architecture> <period> <out_csv_dir> [ks=4,5,6,7,8,9] [sexes=female,male]
a <- commandArgs(trailingOnly = TRUE)
scope <- a[1]; arch <- a[2]; period <- a[3]; OUT <- a[4]
ks <- if (length(a) >= 5) as.integer(strsplit(a[5], ",")[[1]]) else 4:9
sexes <- if (length(a) >= 6) strsplit(a[6], ",")[[1]] else c("female", "male")
EX <- Sys.getenv("V2_EX", file.path(Sys.getenv("HOME"), "mnt/revision/variables/GAMM_day_night_spatial_20261028/saturation_intervals_V2/execution_20261102"))
CODE <- file.path(EX, "code")
source(file.path(CODE, "v2_model_common.R"), encoding = "UTF-8")
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x
dir.create(OUT, showWarnings = FALSE, recursive = TRUE); dir.create(file.path(OUT, "curves"), showWarnings = FALSE)
PROG <- file.path(OUT, "KLOG60_PROGRESS.csv")
logline <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), ..., "\n")
job <<- list(op = "fit", scope = scope, architecture = arch, city = "pooled", period = period, family = "quasipoisson_log",
             k = 6L, fold = "all", rep = 0L, root = EX, variant = "expanded")
for (e in parse(file.path(CODE, "v2_fit.R"), encoding = "UTF-8"))
  if (is.call(e) && identical(e[[1]], as.name("<-")) && identical(as.character(e[[2]]), "prepare_sex")) eval(e, globalenv())
load_wide <- function() {
  cache <- file.path(EX, "cache/raw", paste0(scope, "__", period, "__pooled__variant_expanded"))
  mark <- jsonlite::fromJSON(file.path(cache, "RAW_WIDE.json"), simplifyVector = FALSE)
  sha <- v2_sha(file.path(cache, "RAW_WIDE.rds"))   # streams bytes; avoids NAS owner-SID lookup that breaks digest(file=) on Windows
  if (!identical(sha, mark$cache$sha256)) stop("RAW_WIDE.rds SHA differs from its marker")
  wide <- readRDS(file.path(cache, "RAW_WIDE.rds")); setDT(wide)
  meta <- readRDS(file.path(EX, "inputs/expanded", scope, period, "METADATA.rds"))
  scaling <- as.data.table(meta$scaling)[match(V2_FIELDS, field)]
  att <- v2_attach_folds(wide, EX, FALSE); wide <- att$data
  wide <- v2_scale(wide, scaling); wide[, c(V2_FIELDS, intersect("spatial_block_5km", names(wide))) := NULL]
  wide[, case_weight := v2_weights(wide, 0L, att$blocks)]
  list(wide = wide, scaling = scaling, sha = sha)
}
for (sex in sexes) {
  logline("LOAD", scope, arch, period, sex)
  W <- load_wide(); sc <- W$scaling[field == "road_mean"]
  z_all <- W$wide$z_road_mean
  prep <- prepare_sex(W$wide, sex); W$wide <- NULL; gc()
  if (prep$status != "ready") { fwrite(data.table(scope, arch, period, sex, k = NA, status = prep$status), PROG, append = file.exists(PROG)); next }
  d <- as.data.frame(prep$data); prep$data <- NULL; gc()
  raw <- d$z_road_mean * sc$sample_sd + sc$center
  if (any(!is.finite(raw) | raw <= 0)) stop("non-positive L; log scale undefined")
  d$xL <- log(raw); uL <- sort(unique(log(z_all * sc$sample_sd + sc$center)))   # production rule: unique training L of the pooled frame
  y <- d$y; pw <- d$case_weight
  qz <- quantile(d$z_road_mean, c(0, .005, .01, .99, .995, 1), type = 7)
  for (k in ks) {
    t0 <- proc.time()[3]; w <- character(); err <- NA_character_
    kn <- as.numeric(quantile(uL, seq(0, 1, length.out = k), names = FALSE, type = 7))
    if (Sys.getenv("KNOTS") == "even") kn <- c(seq(min(d$xL), quantile(d$xL, .99), length.out = k - 1), max(d$xL))   # 1110: evenly spaced on log L (min..P99) + boundary at max
    f0 <- v2_formula(arch, k, prep$controls, length(prep$city_levels))
    f <- as.formula(gsub("z_road_mean", "xL", paste(deparse(f0, width.cutoff = 500L), collapse = " "), fixed = TRUE), env = environment(f0))
    m <- tryCatch(withCallingHandlers(bam(f, data = d, weights = case_weight, family = quasipoisson(link = "log"), method = "fREML",
            discrete = 10000L, nthreads = 2L, gc.level = 2L, chunk.size = 5000L, na.action = na.fail, drop.unused.levels = FALSE,
            knots = list(xL = kn), control = gam.control(trace = FALSE, epsilon = 1e-7, maxit = 200L)),
            warning = function(cw) { w <<- unique(c(w, conditionMessage(cw))); invokeRestart("muffleWarning") }),
          error = function(e) { err <<- conditionMessage(e); NULL })
    secs <- unname(proc.time()[3] - t0)
    if (is.null(m)) {
      fwrite(data.table(scope = scope, architecture = arch, period = period, sex = sex, k = k, status = "error", error = err, secs = secs),
             PROG, append = file.exists(PROG)); logline("ERROR", k, err); next }
    mu <- fitted(m); pc <- pw * (y - mu)^2 / mu; edf <- m$edf; sp <- m$sp
    conv <- isTRUE(m$mgcv.conv) && isTRUE(m$converged) && !any(grepl("converge|non-finite|step failure", w))
    finite <- all(is.finite(coef(m))) && all(is.finite(m$Vp)) && all(is.finite(mu))
    phi_p <- sum(pc) / (length(y) - sum(edf))
    r <- data.table(scope = scope, architecture = arch, period = period, sex = sex, k = k, status = if (conv && finite) "accepted" else "rejected",
      converged = conv, finite = finite, iterations = m$iter %||% NA_integer_, secs = round(secs, 1), n = length(y), n_positive = sum(y > 0),
      phi_reported = m$sig2, phi_pearson = phi_p, phi_ratio = m$sig2 / phi_p,
      phi_deviance = sum(quasipoisson()$dev.resids(y, mu, pw)) / (length(y) - sum(edf)),
      top10_X2_share = sum(sort(pc, decreasing = TRUE)[1:10]) / sum(pc), n_pos_mu_lt_1e6 = sum(y > 0 & mu < 1e-6), min_mu_pos = min(mu[y > 0]),
      edf_total = sum(edf), edf_road = sum(edf[grepl("^s\\(xL\\)", names(edf))]), edf_citydev = sum(edf[grepl("^s\\(city_fac", names(edf))]),
      sp_road = unname(sp["s(xL)"]) %||% NA_real_, max_outer_grad = if (length(m$outer.info$grad)) max(abs(m$outer.info$grad)) else NA_real_,
      deviance = m$deviance, null_deviance = m$null.deviance, knots_logL = paste(signif(kn, 6), collapse = ";"),
      warnings = paste(w, collapse = " | "), mgcv = as.character(packageVersion("mgcv")), raw_sha256 = W$sha)
    r[, tail_healthy := n_pos_mu_lt_1e6 == 0 & phi_ratio > 0.5 & phi_ratio < 2]
    fwrite(r, PROG, append = file.exists(PROG))
    ## shared-curve on the original L axis (centred at median L of the frame)
    g <- sort(unique(c(seq(qz[1], qz[6], length.out = 250), seq(qz[2], qz[5], length.out = 250))))
    nd <- d[rep(1L, length(g)), , drop = FALSE]; nd$xL <- log(g * sc$sample_sd + sc$center)
    pp <- predict(m, nd, type = "terms", se.fit = TRUE, discrete = FALSE); j <- which(colnames(pp$fit) == "s(xL)")
    ref <- approx(g, pp$fit[, j], xout = median(d$z_road_mean))$y
    fwrite(data.table(scope = scope, architecture = arch, period = period, sex = sex, k = k, z_road_mean = g, L = g * sc$sample_sd + sc$center,
                      in_P005_P995 = g >= qz[2] & g <= qz[5], in_P01_P99 = g >= qz[3] & g <= qz[4], f = pp$fit[, j] - ref, se = pp$se.fit[, j]),
           file.path(OUT, "curves", sprintf("KLOG_curve_%s_%s_%s_k%d_%s.csv", scope, arch, period, k, sex)))
    if (arch == "city_deviation" && Sys.getenv("CITY_CURVES", "0") == "1") {   # city-specific curves for regime/clustering stability
      cc <- rbindlist(lapply(prep$city_levels, function(cty) {
        zc <- d$z_road_mean[d$city == cty]; if (length(zc) < 50) return(NULL)
        gz <- seq(quantile(zc, .01), quantile(zc, .99), length.out = 60)
        ndc <- d[match(cty, d$city), , drop = FALSE][rep(1L, 60), , drop = FALSE]; ndc$xL <- log(gz * sc$sample_sd + sc$center)
        tt <- predict(m, ndc, type = "terms", discrete = FALSE); jj <- grep("xL", colnames(tt))
        data.table(scope, architecture = arch, period, sex, k, city = cty, L = gz * sc$sample_sd + sc$center, f = rowSums(tt[, jj, drop = FALSE]))
      }))
      fwrite(cc, file.path(OUT, "curves", sprintf("KLOG_citycurves_%s_%s_k%d_%s.csv", arch, period, k, sex)))
    }
    logline("FIT", scope, arch, period, sex, "k", k, r$status, "it", r$iterations, "secs", round(secs), "phi", signif(r$phi_reported, 4), "tail", r$tail_healthy)
    rm(m, mu, pc); gc()
  }
  rm(d, W, prep); gc()
}
logline("GROUP_DONE", scope, arch, period)
