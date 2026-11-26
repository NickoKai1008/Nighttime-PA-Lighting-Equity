## 53_citydev_curves_v2.R — city-specific light curves "same 25-city shape, mild city differences":
## global s(log L) (k=8, KLOG) + penalised city deviations s(city, log L, bs="sz") with a SMALLER basis (default k_dev=5) than the global smooth,
## production controls and city-deviation design otherwise (v2_formula). For every city: curve on its own P2.5–P97.5 light range, anchored at the
## city's median light, with model-based pointwise 95% bands (lpmatrix x Vp over the global + deviation columns), the city level K
## (mean fitted PA per grid-day at median-lit roads, other covariates as observed) and support counts. Also the global (city-average) curve.
## Usage: Rscript 53_citydev_curves_v2.R <scope> <period> <out_dir> [k_dev=5] [sexes]
a <- commandArgs(trailingOnly = TRUE); scope <- a[1]; period <- a[2]; OUT <- a[3]; KDEV <- if (length(a) >= 4) as.integer(a[4]) else 5L
sexes <- if (length(a) >= 5) strsplit(a[5], ",")[[1]] else c("female", "male"); k <- as.integer(a[6]); stopifnot(length(k) == 1L, k %in% c(4L, 5L, 6L, 7L, 9L))
EX <- Sys.getenv("V2_EX"); CODE <- file.path(EX, "code"); source(file.path(CODE, "v2_model_common.R"), encoding = "UTF-8")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE); logline <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), ..., "\n")
job <<- list(op = "fit", scope = scope, architecture = "city_deviation", city = "pooled", period = period, family = "quasipoisson_log", k = k, fold = "all", rep = 0L, root = EX, variant = "expanded")
for (e in parse(file.path(CODE, "v2_fit.R"), encoding = "UTF-8"))
  if (is.call(e) && identical(e[[1]], as.name("<-")) && identical(as.character(e[[2]]), "prepare_sex")) eval(e, globalenv())
cache <- file.path(EX, "cache/raw", paste0(scope, "__", period, "__pooled__variant_expanded"))
mark <- jsonlite::fromJSON(file.path(cache, "RAW_WIDE.json"), simplifyVector = FALSE); if (!identical(v2_sha(file.path(cache, "RAW_WIDE.rds")), mark$cache$sha256)) stop("RAW_WIDE SHA mismatch")
for (sex in sexes) {
  wide <- readRDS(file.path(cache, "RAW_WIDE.rds")); setDT(wide)
  TC <- Sys.getenv("TEST_CITIES", ""); if (nzchar(TC)) wide <- wide[city %in% strsplit(TC, ",")[[1]]]   # for quick tests only
  logline("LOADED", nrow(wide))
  meta <- readRDS(file.path(EX, "inputs/expanded", scope, period, "METADATA.rds")); sc <- as.data.table(meta$scaling)[match(V2_FIELDS, field)]
  att <- v2_attach_folds(wide, EX, FALSE); wide <- att$data; Lraw <- wide$road_mean
  wide <- v2_scale(wide, sc); wide[, c(V2_FIELDS, intersect("spatial_block_5km", names(wide))) := NULL]; wide[, z_road_mean := log(Lraw)]; wide[, case_weight := 1]
  uL <- sort(unique(wide$z_road_mean)); prep <- prepare_sex(wide, sex); rm(wide, att, Lraw, meta, sc); gc(); d <- as.data.frame(prep$data); prep$data <- NULL; gc()
  XS <- Sys.getenv("XSCALE", "log"); trf <- switch(XS, log = log, sqrt = sqrt, lin = function(L) L)   # 1110: basis scale for the light smooths
  d$zlog <- d$z_road_mean; d$xL <- trf(exp(d$zlog)); d$xLd <- d$xL
  kn <- as.numeric(quantile(uL, seq(0, 1, length.out = k), names = FALSE, type = 7)); knd <- as.numeric(quantile(uL, seq(0, 1, length.out = KDEV), names = FALSE, type = 7))
  if (Sys.getenv("KNOTS") == "even") { zq <- quantile(d$xL, .99); kn <- c(seq(min(d$xL), zq, length.out = k - 1), max(d$xL)); knd <- c(seq(min(d$xL), zq, length.out = KDEV - 1), max(d$xL)) }   # 1110 option
  f0 <- paste(deparse(v2_formula("city_deviation", k, prep$controls, length(prep$city_levels)), width.cutoff = 500L), collapse = " ")
  f0 <- sub(sprintf("s(city_fac, z_road_mean, k = %d", k), sprintf("s(city_fac, xLd, k = %d", KDEV), f0, fixed = TRUE)
  f0 <- gsub("z_road_mean", "xL", f0, fixed = TRUE); if (!grepl("xLd", f0)) stop("deviation term not found in formula: ", f0)
  t0 <- proc.time()[3]; wn <- character()
  m <- withCallingHandlers(bam(as.formula(f0, env = asNamespace("mgcv")), data = d, weights = case_weight, family = quasipoisson(link = "log"), method = "fREML",
         discrete = 10000L, nthreads = 2L, gc.level = 2L, chunk.size = 5000L, na.action = na.fail, drop.unused.levels = FALSE,
         knots = list(xL = kn, xLd = knd), control = gam.control(epsilon = 1e-7, maxit = 200L)), warning = function(cw) { wn <<- unique(c(wn, conditionMessage(cw))); invokeRestart("muffleWarning") })
  conv <- isTRUE(m$mgcv.conv) && isTRUE(m$converged) && !any(grepl("converge|non-finite", wn)); mu <- fitted(m)
  logline("FIT", scope, period, sex, "k_dev", KDEV, if (conv) "accepted" else "REJECTED", "it", m$iter, "secs", round(proc.time()[3] - t0))
  b <- coef(m); V <- m$Vp; cn <- names(b); cg <- grep("^s\\(xL\\)", cn); cd <- grep("^s\\(city_fac,xLd\\)", cn)
  curve <- function(rowdat, Lg, Lref, cols) {
    nd <- rowdat[rep(1L, length(Lg) + 1L), , drop = FALSE]; nd$xL <- nd$xLd <- trf(c(Lg, Lref))
    X <- predict(m, nd, type = "lpmatrix", discrete = FALSE)[, cols, drop = FALSE]; D <- sweep(X[seq_along(Lg), , drop = FALSE], 2, X[length(Lg) + 1L, ])
    list(f = as.numeric(D %*% b[cols]), se = sqrt(pmax(rowSums((D %*% V[cols, cols]) * D), 0)))
  }
  out <- list(); bc <- list()
  for (cty in prep$city_levels) {
    ii <- which(d$city == cty); if (length(ii) < 100) next
    Lc <- exp(d$zlog[ii]); Lg <- seq(quantile(Lc, .025), quantile(Lc, .975), length.out = 80); Lm <- median(Lc)
    cv <- curve(d[ii[1], , drop = FALSE], Lg, Lm, c(cg, cd))
    fi <- approx(log(Lg), cv$f, xout = log(Lc), rule = 2)$y; K <- mean(mu[ii] * exp(-fi))   # PA per grid-day at the city's median light
    out[[cty]] <- data.table(scope, period, sex, k, k_dev = KDEV, city = cty, L = Lg, f = cv$f, se = cv$se, lo = cv$f - 1.96 * cv$se, hi = cv$f + 1.96 * cv$se,
                             L_median = Lm, K_PA_per_griddday_at_median = K, n_rows = length(ii), n_positive = sum(d$y[ii] > 0), converged = conv)
    if (Sys.getenv("BINCHECK") == "1") {   # 1111: per-grid sums for the fit check - y, fitted mu, and m0 = mu with the city's light effect set to its median light
      Lgf <- exp(seq(log(min(Lc)), log(max(Lc)), length.out = 200)); cvf <- curve(d[ii[1], , drop = FALSE], Lgf, Lm, c(cg, cd))
      gi <- approx(log(Lgf), cvf$f, xout = log(Lc), rule = 2)$y
      bc[[cty]] <- data.table(city = cty, grid_id = if ("grid_id" %in% names(d)) as.character(d$grid_id[ii]) else NA_character_, logL = d$zlog[ii], y = d$y[ii], mu = mu[ii], m0 = mu[ii] * exp(-gi))[
        , .(n = .N, y = sum(y), mu = sum(mu), m0 = sum(m0), logL = mean(logL)), by = .(city, grid_id)]
    }
  }
  fwrite(rbindlist(out), file.path(OUT, sprintf("CITY_CURVES_V2_%s_%s_kdev%d_%s.csv", scope, period, KDEV, sex)))
  if (length(bc)) fwrite(rbindlist(bc), file.path(OUT, sprintf("BINCHECK_GRID_%s_%s_kdev%d_%s.csv.gz", scope, period, KDEV, sex)))
  Lall <- exp(d$zlog); Lg <- seq(quantile(Lall, .01), quantile(Lall, .99), length.out = 200); cv <- curve(d[1, , drop = FALSE], Lg, median(Lall), cg)
  fwrite(data.table(scope, period, sex, k, k_dev = KDEV, L = Lg, f = cv$f, se = cv$se, L_median = median(Lall), K_PA_per_griddday_at_median = mean(mu * exp(-approx(log(Lg), cv$f, xout = log(Lall), rule = 2)$y))),
         file.path(OUT, sprintf("GLOBAL_CURVE_V2_%s_%s_kdev%d_%s.csv", scope, period, KDEV, sex)))
  fwrite(data.table(scope, period, sex, k_dev = KDEV, status = if (conv) "accepted" else "rejected", iterations = m$iter, deviance = m$deviance, edf_total = sum(m$edf),
                    edf_dev = sum(m$edf[grepl("city_fac,xLd", names(m$edf))]), warnings = paste(wn, collapse = " | ")), file.path(OUT, "PROGRESS.csv"), append = file.exists(file.path(OUT, "PROGRESS.csv")))
  rm(m, d, mu, prep); gc()
}
