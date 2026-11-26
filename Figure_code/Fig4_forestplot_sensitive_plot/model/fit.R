# Refit the two all-city nighttime models used in panel a.
library(mgcv)

script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
root <- dirname(normalizePath(script))
input <- file.path(root, "data")
fig_dir <- dirname(root)

observations <- read.csv(file.path(input, "observations.csv.gz"),
                         colClasses = c("factor", "integer", "integer", "numeric", "numeric"))
cells <- read.csv(file.path(input, "cells.csv.gz"))
periods <- read.csv(file.path(input, "periods.csv.gz"))
cities <- read.csv(file.path(input, "cities.csv"))
scaling <- read.csv(file.path(input, "scaling.csv"))
reference <- read.csv(file.path(fig_dir, "data", "coefficients.csv"))

# The within-city keys link each observation to its static cell and daily predictors.
index_within_city <- function(city, key, lookup_city, lookup_key) {
  result <- integer(length(key))
  for (name in unique(city)) {
    rows <- which(city == name)
    lookup <- which(lookup_city == name)
    result[rows] <- lookup[match(key[rows], lookup_key[lookup])]
  }
  stopifnot(!anyNA(result))
  result
}
cell_id <- index_within_city(observations$city, observations$cell, cells$city, cells$cell)
period_id <- index_within_city(observations$city, observations$period, periods$city, periods$period)
city_id <- match(observations$city, cities$city)

d <- data.frame(city_fac = factor(observations$city),
                date_num = periods$date_num[period_id],
                is_weekend = periods$is_weekend[period_id],
                z_road_mean = log(cells$road_mean[cell_id]))
for (i in seq_len(nrow(scaling))) {
  field <- scaling$field[i]
  values <- if (field %in% names(cells)) cells[[field]][cell_id] else periods[[field]][period_id]
  d[[paste0("z_", field)]] <- (values - scaling$center[i]) / scaling$sample_sd[i]
}
city_terms <- setdiff(names(cities), "city")
for (field in city_terms) d[[paste0("c_", field)]] <- as.numeric(scale(cities[[field]]))[city_id]

weather <- paste0("z_", c("tavg", "prcp", "wspd", "pres", "rhum"))
grid <- setdiff(paste0("z_", scaling$field), c("z_road_relative_iqr", weather))
rhs <- c(grid, weather, "is_weekend", "s(date_num,k=5,bs='cr')",
         "z_road_relative_iqr", paste0("c_", city_terms),
         "z_road_mean", "s(city_fac,bs='re')")
formula <- as.formula(paste("y ~", paste(rhs, collapse = " + ")))
rm(values, cells, periods, cell_id, period_id, city_id)
stopifnot(nrow(d) == 2800643L, nlevels(d$city_fac) == 25L)

refit <- function(sex) {
  message("Fitting ", sex, " nighttime activity")
  d$y <- observations[[sex]]
  fit <- bam(formula, data = d, family = quasipoisson(link = "log"), method = "fREML",
             discrete = 10000L, nthreads = 2L, gc.level = 2L, chunk.size = 5000L,
             na.action = na.fail, drop.unused.levels = FALSE,
             control = gam.control(epsilon = 1e-7, maxit = 200L))
  stopifnot(isTRUE(fit$converged))
  out <- reference[reference$sex == sex, ]
  ptab <- summary(fit)$p.table[out$term, , drop = FALSE]
  out$estimate <- ptab[, 1]
  out$se <- ptab[, 2]
  out$p <- ptab[, 4]
  out$lo_95_normal <- out$estimate - qnorm(0.975) * out$se
  out$hi_95_normal <- out$estimate + qnorm(0.975) * out$se
  out
}
fitted <- do.call(rbind, lapply(c("female", "male"), refit))
output <- file.path(fig_dir, "output")
dir.create(output, showWarnings = FALSE)
write.csv(fitted, file.path(output, "refit_coefficients.csv"), row.names = FALSE)
for (field in c("estimate", "se", "p", "lo_95_normal", "hi_95_normal")) {
  original <- reference[[field]][match(paste(fitted$sex, fitted$term), paste(reference$sex, reference$term))]
  delta <- max(abs(fitted[[field]] - original))
  message(field, " — maximum absolute difference: ", signif(delta, 4))
  if (!is.finite(delta) || delta > 1e-4) stop("Check the refit against the panel-a source table.")
}
