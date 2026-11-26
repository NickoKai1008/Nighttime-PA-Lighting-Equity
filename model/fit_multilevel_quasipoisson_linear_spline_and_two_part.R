## 27_multilevel_main.R — pre-specified (2026-10-21 six-indicator framework) cross-city multilevel specification.
## City random intercept per sex-specific fit (s(city_fac, bs="re")) instead of city fixed effects, so the 12 city-level
## constants are identifiable. L enters as s(log L) (KLOG). age35plus_reported_share is NOT used (activity-derived).
## Specs (A5 nested sequence + FE reference):
##   B        : 11 grid controls (excl. H) + 5 weather + is_weekend + s(date) + 8 social city-level + city RE   [no lighting]
##   B_L      : B + s(log L)
##   B_L_H    : B + s(log L) + H
##   B_L_S    : B + s(log L) + alpha, D, c, ULPI
##   FULL     : B + s(log L) + H + alpha, D, c, ULPI        (main specification)
##   FE_REF   : production FE design (0 + city_fac) + H + 11 grid + weather + s(log L)   (city constants absorbed)
## Usage: Rscript 27_multilevel_main.R <scope> <period> <k> <specs comma> <out_dir> <city_level_csv> [sexes]
a <- commandArgs(trailingOnly = TRUE)
scope <- a[1]; period <- a[2]; k <- as.integer(a[3]); specs <- strsplit(a[4], ",")[[1]]; OUT <- a[5]; CL <- a[6]
sexes <- if (length(a) >= 7) strsplit(a[7], ",")[[1]] else c("female", "male")
EX <- Sys.getenv("V2_EX", file.path(Sys.getenv("HOME"), "mnt/revision/variables/GAMM_day_night_spatial_20261028/saturation_intervals_V2/execution_20261102"))
CODE <- file.path(EX, "code"); source(file.path(CODE, "v2_model_common.R"), encoding = "UTF-8")
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x
dir.create(file.path(OUT, "curves"), showWarnings = FALSE, recursive = TRUE); PR <- file.path(OUT, "MULTILEVEL_PROGRESS.csv"); CO <- file.path(OUT, "MULTILEVEL_CITY_COEF.csv")
logline <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), ..., "\n")
job <<- list(op = "fit", scope = scope, architecture = "shared", city = "pooled", period = period, family = "quasipoisson_log",
             k = k, fold = "all", rep = 0L, root = EX, variant = "expanded")
for (e in parse(file.path(CODE, "v2_fit.R"), encoding = "UTF-8"))
  if (is.call(e) && identical(e[[1]], as.name("<-")) && identical(as.character(e[[2]]), "prepare_sex")) eval(e, globalenv())
STRUCT <- c("ntl_alpha", "ntl_D", "ntl_c", "ulpi")
SOCIAL <- c("female_lfpr", "female_upper_secondary", "wbl_supportive_childcare", "public_open_space_share", "numbeo_safety_walk_night",
            "maternal_mortality_ratio", "adolescent_birth_rate", "female_parliament_share")
cl <- fread(CL); for (v in c(STRUCT, SOCIAL)) cl[, (paste0("c_", v)) := as.numeric(scale(get(v)))]   # z across the 25 cities
cache <- file.path(EX, "cache/raw", paste0(scope, "__", period, "__pooled__variant_expanded"))
mark <- jsonlite::fromJSON(file.path(cache, "RAW_WIDE.json"), simplifyVector = FALSE)
if (!identical(v2_sha(file.path(cache, "RAW_WIDE.rds")), mark$cache$sha256)) stop("RAW_WIDE SHA mismatch")
wide <- readRDS(file.path(cache, "RAW_WIDE.rds")); setDT(wide)
DROP <- Sys.getenv("ML_DROP_CITIES", ""); if (nzchar(DROP)) { DROP <- strsplit(DROP, ",")[[1]]; if (!all(DROP %in% wide$city)) stop("unknown city in ML_DROP_CITIES"); wide <- wide[!city %in% DROP]; logline("DROPPED", paste(DROP, collapse = ",")) }   # influence checks (LOCO/LORO)
SOLF <- Sys.getenv("ML_SOLAR_CSV", ""); if (nzchar(SOLF)) {   # 1110: solar-period outcomes (45_build_solar_outcomes.py; same allocation as the clock replication), covered grids only
  SV <- Sys.getenv("ML_SOLAR_VAR", "solar"); so <- if (grepl("[.]gz$", SOLF) && nzchar(Sys.which("gzip"))) fread(cmd = paste("gzip -dc", shQuote(SOLF)), colClasses = c(grid_id = "character", date = "character")) else fread(SOLF, colClasses = c(grid_id = "character", date = "character"))
  so <- so[, .(city, gid = grid_id, dnum = as.numeric(date_num), yf = get(sprintf("female_%s_%s", SV, period)), ym = get(sprintf("male_%s_%s", SV, period)))]
  wide[, `:=`(gid = as.character(grid_id), dnum = as.numeric(date_num))]; n0 <- nrow(wide)
  wide <- merge(wide, so, by = c("city", "gid", "dnum")); wide[, `:=`(y_female = yf, y_male = ym)]; wide[, c("gid", "dnum", "yf", "ym") := NULL]
  logline("SOLAR outcomes", SV, period, "rows", nrow(wide), "of", n0) }
meta <- readRDS(file.path(EX, "inputs/expanded", scope, period, "METADATA.rds")); sc <- as.data.table(meta$scaling)[match(V2_FIELDS, field)]
att <- v2_attach_folds(wide, EX, FALSE); wide <- att$data; Lraw <- wide$road_mean
wide <- v2_scale(wide, sc); wide[, c(V2_FIELDS, intersect("spatial_block_5km", names(wide))) := NULL]; wide[, z_road_mean := log(Lraw)]
wide[, case_weight := 1]
REPS <- { r <- Sys.getenv("ML_FRW_REPS", ""); if (nzchar(r)) { b <- as.integer(strsplit(r, "-")[[1]]); seq(b[1], b[length(b)]) } else 0L }   # 1109: FRW replicates in-process (2-km block Gamma(1,1) weights, same seed rule as production)
kn <- as.numeric(quantile(sort(unique(wide$z_road_mean)), seq(0, 1, length.out = k), names = FALSE, type = 7))
qz <- quantile(wide$z_road_mean, c(.01, .99)); g <- seq(qz[1], qz[2], length.out = 200)
for (FRW_REP in REPS) {
RSUF <- if (FRW_REP > 0) sprintf("_r%03d", FRW_REP) else ""
if (FRW_REP > 0) { wide[, case_weight := v2_weights(wide, FRW_REP, att$blocks)]; PR <- file.path(OUT, sprintf("MULTILEVEL_PROGRESS%s.csv", RSUF)); CO <- file.path(OUT, sprintf("MULTILEVEL_CITY_COEF%s.csv", RSUF)); logline("FRW rep", FRW_REP) }
for (sex in sexes) {
  prep <- prepare_sex(wide, sex); d <- as.data.frame(prep$data)
  d <- merge(d, as.data.frame(cl[, c("city", paste0("c_", c(STRUCT, SOCIAL))), with = FALSE]), by = "city", all.x = TRUE, sort = FALSE)
  if (anyNA(d[, paste0("c_", c(STRUCT, SOCIAL))])) stop("city-level values missing for some city")
  PPF <- Sys.getenv("ML_PP_CSV", "")   # 1109: PerceptAtlas street-view perception (daytime imagery), NA unless the 200-m grid is fully covered
  if (nzchar(PPF)) { if (!exists("ppd")) ppd <- fread(PPF, colClasses = c(grid_id = "character"))[, .(city, grid_id, p_safer, p_livelier)]
    d$grid_id <- as.character(d$grid_id); d <- merge(d, as.data.frame(ppd), by = c("city", "grid_id"), all.x = TRUE, sort = FALSE)
    cv <- !is.na(d$p_safer); for (v in c("p_safer", "p_livelier")) d[[paste0("z_", v)]] <- (d[[v]] - mean(d[[v]][cv])) / sd(d[[v]][cv])
    logline("PP coverage", sex, round(mean(cv), 3)) }
  grid11 <- setdiff(prep$controls, c("z_road_relative_iqr", V2_WEATHER, paste0("z_", V2_WEATHER), "is_weekend"))
  base_rhs <- c(grid11, paste0("z_", V2_WEATHER)[paste0("z_", V2_WEATHER) %in% names(d)], "is_weekend", "s(date_num,k=5,bs='cr')")
  ## CPTED blocks from EXISTING variables only (no new data). FULL = CP3 + entrapment block.
  CP_DEST <- intersect(c("z_poi_density_500m", "z_poi_type_mix_500m", "z_walkability_global_z"), names(d))
  CP_PROS <- intersect(c("z_isovist_fraction", "z_road_relative_iqr"), names(d))
  CP_CONC <- intersect(c("z_BD200"), names(d))
  CP_GREEN <- intersect(c("z_NDVI"), names(d))   # green exposure kept as its own block (not CPTED concealment)
  CP_ENTR <- intersect(c("z_distance_to_branch_m", "z_mapped_terminal_branch", "z_ssx_choice_r5000"), names(d))
  cp_core <- c(setdiff(base_rhs, c(CP_DEST, CP_PROS, CP_CONC, CP_ENTR, CP_GREEN)))
  if (!file.exists(file.path(OUT, sprintf("GRID_CONTROL_COR_%s_%s_%s.csv", scope, period, sex)))) {
    gc12 <- intersect(c(CP_DEST, CP_PROS, CP_CONC, CP_ENTR, CP_GREEN, "z_Pop_Female_20_50", "z_gdp_2022"), names(d))
    cm <- cor(d[, gc12]); vif <- diag(solve(cm))
    fwrite(data.table(variable = gc12, block = ifelse(gc12 %in% CP_DEST, "destinations", ifelse(gc12 %in% CP_PROS, "prospect", ifelse(gc12 %in% CP_CONC, "concealment", ifelse(gc12 %in% CP_ENTR, "entrapment", ifelse(gc12 %in% CP_GREEN, "green exposure", "population/economy"))))),
      VIF = vif, max_abs_r = apply(abs(cm - diag(nrow(cm))), 1, max), cm), file.path(OUT, sprintf("GRID_CONTROL_COR_%s_%s_%s.csv", scope, period, sex)))
  }
  for (spec in specs) {
    sL <- sprintf("s(z_road_mean,k=%d,bs='cr')", k); sLby <- sprintf("s(z_road_mean,k=%d,bs='cr',by=byf)", k)
    full_rhs <- c("1", base_rhs, "z_road_relative_iqr", paste0("c_", c(STRUCT, SOCIAL)), sL, "s(city_fac,bs='re')")
    rhs <- switch(spec,
      FULL_EVEN  = full_rhs,   # 1110: same as FULL, knots evenly spaced on log L from the minimum to P99 (+ boundary knot at the maximum)
      FULL_SQRT  = full_rhs, FULL_LINL = full_rhs, FULL_LINLQ = full_rhs,   # 1110: light smooth on sqrt(L) / linear L (as the 10-28 atlas) instead of log L
      FULL_DET   = c(setdiff(full_rhs, sL), sprintf("s(z_road_mean,k=%d,bs='cr',by=det)", k)),   # 1110: detection split - below-detection (DN-floor) roads are the reference level; smooth only over detected light
      FULL_PPSUB = full_rhs,                                              # FULL on the perception-covered subsample (baseline for FULL_PP)
      FULL_PP    = c(full_rhs, "z_p_safer"),                              # + perceived safety (daytime street view)
      FULL_PP2   = c(full_rhs, "z_p_safer", "z_p_livelier"),              # + perceived liveliness
      BY_SAFE    = c(setdiff(full_rhs, sL), "byf", sLby),                 # light curve by within-city tercile of perceived safety
      BY_WKND    = c(setdiff(full_rhs, sL), sLby),                        # light curve weekday vs weekend (is_weekend already in base)
      BY_HALF    = c(setdiff(full_rhs, sL), "byf", sLby),                 # light curve first vs second half of the study period
      B     = c("1", base_rhs, paste0("c_", SOCIAL), "s(city_fac,bs='re')"),
      B_L   = c("1", base_rhs, paste0("c_", SOCIAL), sprintf("s(z_road_mean,k=%d,bs='cr')", k), "s(city_fac,bs='re')"),
      B_L_H = c("1", base_rhs, "z_road_relative_iqr", paste0("c_", SOCIAL), sprintf("s(z_road_mean,k=%d,bs='cr')", k), "s(city_fac,bs='re')"),
      B_L_S = c("1", base_rhs, paste0("c_", c(STRUCT, SOCIAL)), sprintf("s(z_road_mean,k=%d,bs='cr')", k), "s(city_fac,bs='re')"),
      FULL  = c("1", base_rhs, "z_road_relative_iqr", paste0("c_", c(STRUCT, SOCIAL)), sprintf("s(z_road_mean,k=%d,bs='cr')", k), "s(city_fac,bs='re')"),
      FE_REF = c("0+city_fac", base_rhs, "z_road_relative_iqr", sprintf("s(z_road_mean,k=%d,bs='cr')", k)),
      FULL_noMMR    = c("1", base_rhs, "z_road_relative_iqr", paste0("c_", c(STRUCT, setdiff(SOCIAL, "maternal_mortality_ratio"))), sprintf("s(z_road_mean,k=%d,bs='cr')", k), "s(city_fac,bs='re')"),
      FULL_noNUMBEO = c("1", base_rhs, "z_road_relative_iqr", paste0("c_", c(STRUCT, setdiff(SOCIAL, "numbeo_safety_walk_night"))), sprintf("s(z_road_mean,k=%d,bs='cr')", k), "s(city_fac,bs='re')"),
      CP0 = c("1", cp_core, paste0("c_", c(STRUCT, SOCIAL)), sprintf("s(z_road_mean,k=%d,bs='cr')", k), "s(city_fac,bs='re')"),
      CP1 = c("1", cp_core, CP_DEST, paste0("c_", c(STRUCT, SOCIAL)), sprintf("s(z_road_mean,k=%d,bs='cr')", k), "s(city_fac,bs='re')"),
      CP2 = c("1", cp_core, CP_DEST, CP_PROS, paste0("c_", c(STRUCT, SOCIAL)), sprintf("s(z_road_mean,k=%d,bs='cr')", k), "s(city_fac,bs='re')"),
      CP3 = c("1", cp_core, CP_DEST, CP_PROS, CP_CONC, paste0("c_", c(STRUCT, SOCIAL)), sprintf("s(z_road_mean,k=%d,bs='cr')", k), "s(city_fac,bs='re')"),
      CP4 = c("1", cp_core, CP_DEST, CP_PROS, CP_CONC, CP_ENTR, paste0("c_", c(STRUCT, SOCIAL)), sprintf("s(z_road_mean,k=%d,bs='cr')", k), "s(city_fac,bs='re')"),
      FULL_LIN = c("1", base_rhs, "z_road_relative_iqr", paste0("c_", c(STRUCT, SOCIAL)), "z_road_mean", "s(city_fac,bs='re')"),   # log L as one linear term (LMM-like single coefficient)
      TP_OCC = c("1", base_rhs, "z_road_relative_iqr", paste0("c_", c(STRUCT, SOCIAL)), sprintf("s(z_road_mean,k=%d,bs='cr')", k), "s(city_fac,bs='re')"),
      TP_POS = c("1", base_rhs, "z_road_relative_iqr", paste0("c_", c(STRUCT, SOCIAL)), sprintf("s(z_road_mean,k=%d,bs='cr')", k), "s(city_fac,bs='re')"))
    fm <- as.formula(paste("y ~", paste(rhs, collapse = " + ")), env = asNamespace("mgcv"))
    t0 <- proc.time()[3]; wn <- character(); fam <- quasipoisson(link = "log"); d0 <- d
    if (spec %in% c("FULL_PPSUB", "FULL_PP", "FULL_PP2", "BY_SAFE")) d <- d0[!is.na(d0$p_safer), , drop = FALSE]
    if (spec == "BY_SAFE") { r <- ave(d$p_safer, d$city, FUN = function(x) rank(x, ties.method = "average") / length(x)); d$byf <- cut(r, c(0, 1/3, 2/3, 1), labels = c("low", "mid", "high"), include.lowest = TRUE) }
    if (spec == "FULL_EVEN") { zq <- quantile(d$z_road_mean, .99); kn_even <- c(seq(min(d$z_road_mean), zq, length.out = k - 1), max(d$z_road_mean)); logline("EVEN knots L", paste(signif(exp(kn_even), 3), collapse = ";")) }
    TRS <- c("FULL_SQRT", "FULL_LINL", "FULL_LINLQ")
    if (spec %in% TRS) { trf <- if (spec == "FULL_SQRT") function(L) sqrt(L) else function(L) L; d$z_log <- d$z_road_mean; d$z_road_mean <- trf(exp(d$z_log))
      kn_t <- if (spec == "FULL_LINLQ") as.numeric(quantile(sort(unique(d$z_road_mean)), seq(0, 1, length.out = k), names = FALSE, type = 7)) else c(seq(min(d$z_road_mean), quantile(d$z_road_mean, .99), length.out = k - 1), max(d$z_road_mean))
      logline("XSCALE knots L", spec, paste(signif(if (spec == "FULL_SQRT") kn_t^2 else kn_t, 3), collapse = ";")) }
    if (spec == "FULL_DET") { Lr <- exp(d$z_road_mean); d$det <- as.numeric(Lr > 1.01 * min(Lr)); kn_det <- as.numeric(quantile(sort(unique(d$z_road_mean[d$det == 1])), seq(0, 1, length.out = k), names = FALSE, type = 7)) }
    if (spec == "BY_WKND") d$byf <- factor(ifelse(d$is_weekend > 0.5, "weekend", "weekday"), levels = c("weekday", "weekend"))
    if (spec == "BY_HALF") d$byf <- factor(ifelse(d$date_num <= median(d$date_num), "early", "late"), levels = c("early", "late"))
    if (spec == "TP_OCC") { d <- d0; d$y <- as.numeric(d0$y > 0); fam <- binomial(link = "logit") }      # two-part: occurrence
    if (spec == "TP_POS") { d <- d0[d0$y > 0, , drop = FALSE]; d$y <- log(d$y); fam <- gaussian() }        # two-part: log positive PA
    m <- withCallingHandlers(bam(fm, data = d, weights = case_weight, family = fam, method = "fREML", discrete = 10000L,
           nthreads = 2L, gc.level = 2L, chunk.size = 5000L, na.action = na.fail, drop.unused.levels = FALSE,
           knots = if (spec %in% TRS) list(z_road_mean = kn_t) else if (spec == "FULL_EVEN") list(z_road_mean = kn_even) else if (spec == "FULL_DET") list(z_road_mean = kn_det) else if (spec != "FULL_LIN" && (grepl("L", spec) || spec %in% c("FULL", "FE_REF") || grepl("^BY_", spec))) list(z_road_mean = kn) else NULL,
           control = gam.control(epsilon = 1e-7, maxit = as.integer(Sys.getenv("ML_MAXIT", "200")))), warning = function(cw) { wn <<- unique(c(wn, conditionMessage(cw))); invokeRestart("muffleWarning") })
    conv <- isTRUE(m$mgcv.conv) && isTRUE(m$converged) && !any(grepl("converge|non-finite", wn))
    mu <- fitted(m); pc <- if (spec %in% c("TP_OCC", "TP_POS")) rep(NA_real_, length(mu)) else (d$y - mu)^2 / mu
    fwrite(data.table(scope, period, sex, k, spec, status = if (conv) "accepted" else "rejected", iterations = m$iter %||% NA, secs = round(unname(proc.time()[3] - t0), 1),
      n = nrow(d), deviance = m$deviance, null_deviance = m$null.deviance, phi = m$sig2, phi_ratio = m$sig2 / (sum(pc) / (nrow(d) - sum(m$edf))),
      edf_total = sum(m$edf), edf_road = sum(m$edf[grepl("z_road_mean", names(m$edf))]), edf_city_re = sum(m$edf[grepl("city_fac", names(m$edf))]),
      n_pos_mu_lt_1e6 = if (spec %in% c("TP_OCC", "TP_POS")) NA_integer_ else sum(d$y > 0 & mu < 1e-6), warnings = paste(wn, collapse = " | ")), PR, append = file.exists(PR))
    cf <- summary(m)$p.table; cf <- data.table(term = rownames(cf), cf)[grepl("^c_|^z_", term)]   # 1109: all parametric controls (incl. z_NDVI), not only city-level + H
    if (nrow(cf)) { setnames(cf, c("term", "estimate", "se", "t", "p")); cf[, `:=`(scope = scope, period = period, sex = sex, k = k, spec = spec)]; fwrite(cf, CO, append = file.exists(CO)) }
    if (spec == "FULL_DET") {
      gd <- seq(min(d$z_road_mean[d$det == 1]), qz[2], length.out = 200); nd <- d[rep(1L, length(gd)), , drop = FALSE]; nd$z_road_mean <- gd; nd$det <- 1
      pp <- predict(m, nd, type = "terms", se.fit = TRUE, discrete = FALSE); j <- which(colnames(pp$fit) == "s(z_road_mean):det")
      fwrite(data.table(scope, period, sex, k, spec, floor_L = min(exp(d$z_road_mean)), share_rows_floor = mean(d$det == 0), L = exp(gd), f = pp$fit[, j], se = pp$se.fit[, j]),
             file.path(OUT, "curves", sprintf("ML_%s_%s_k%d_%s_%s%s.csv", scope, period, k, sex, spec, RSUF)))   # f = log PA contrast vs below-detection roads (floor = 0)
    } else if (grepl("^BY_", spec)) {
      for (lv in levels(d$byf)) { nd <- d[rep(1L, length(g)), , drop = FALSE]; nd$z_road_mean <- g; nd$byf <- factor(lv, levels = levels(d$byf))
        pp <- predict(m, nd, type = "terms", se.fit = TRUE, discrete = FALSE); j <- which(colnames(pp$fit) == paste0("s(z_road_mean):byf", lv))
        ref <- approx(g, pp$fit[, j], xout = median(d$z_road_mean[d$byf == lv]))$y
        fwrite(data.table(scope, period, sex, k, spec, level = lv, n_rows = sum(d$byf == lv), L = exp(g), f = pp$fit[, j] - ref, se = pp$se.fit[, j]),
               file.path(OUT, "curves", sprintf("ML_%s_%s_k%d_%s_%s_%s%s.csv", scope, period, k, sex, spec, lv, RSUF))) }
    } else if (any(grepl("s(z_road_mean", names(m$edf), fixed = TRUE))) {
      nd <- d[rep(1L, length(g)), , drop = FALSE]; nd$z_road_mean <- if (spec %in% TRS) trf(exp(g)) else g
      pp <- predict(m, nd, type = "terms", se.fit = TRUE, discrete = FALSE); j <- which(colnames(pp$fit) == "s(z_road_mean)")
      ref <- approx(g, pp$fit[, j], xout = if (spec %in% TRS) median(d$z_log) else median(d$z_road_mean))$y
      fwrite(data.table(scope, period, sex, k, spec, L = exp(g), f = pp$fit[, j] - ref, se = pp$se.fit[, j]), file.path(OUT, "curves", sprintf("ML_%s_%s_k%d_%s_%s%s.csv", scope, period, k, sex, spec, RSUF)))
    }
    if (spec %in% c("FULL") && Sys.getenv("ML_GRIDAGG", "1") == "1") {
      dir.create(file.path(OUT, "grid_agg"), showWarnings = FALSE)
      dd <- data.table(city = d$city, grid_id = as.character(d$grid_id), y = d$y, mu = mu, pr = (d$y - mu) / sqrt(m$sig2 * mu), logL = d$z_road_mean)
      stat <- intersect(c(CP_DEST, CP_PROS, CP_CONC, CP_ENTR, CP_GREEN, "z_Pop_Female_20_50", "z_gdp_2022"), names(d)); for (v in stat) dd[[v]] <- d[[v]]
      for (v in paste0("z_", V2_WEATHER)) if (v %in% names(d)) dd[[v]] <- d[[v]]
      ag <- dd[, c(list(n_days = .N, y_mean = mean(y), y_pos_share = mean(y > 0), mu_mean = mean(mu), pr_mean = mean(pr)),
                   lapply(.SD, mean)), by = .(city, grid_id), .SDcols = c("logL", stat, intersect(paste0("z_", V2_WEATHER), names(dd)))]
      fwrite(ag, file.path(OUT, "grid_agg", sprintf("GRIDAGG_%s_%s_%s_FULL.csv.gz", scope, period, sex)))
    }
    logline("ML", scope, period, sex, spec, if (conv) "accepted" else "REJECTED", "it", m$iter, "secs", round(unname(proc.time()[3] - t0)), "dev", signif(m$deviance, 7))
    d <- d0; rm(m); gc()
  }
}
}   # FRW replicate loop
