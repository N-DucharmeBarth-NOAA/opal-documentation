#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(RTMB)
  library(dplyr)
  library(tibble)
})

msg <- function(...) {
  cat(format(Sys.time(), "%H:%M:%S"), ..., "\n")
  flush(stdout())
}

max_kkt_gradient <- function(par, grad, lower, upper, tol = 1e-8) {
  active_lower <- is.finite(lower) & abs(par - lower) < tol
  active_upper <- is.finite(upper) & abs(par - upper) < tol
  projected <- grad
  projected[active_lower & projected > 0] <- 0
  projected[active_upper & projected < 0] <- 0
  max(abs(projected), na.rm = TRUE)
}

component_sum <- function(report, name) {
  if (is.null(report[[name]])) return(NA_real_)
  sum(as.numeric(report[[name]]), na.rm = TRUE)
}

sdnr_summary <- function(x, conf = 0.95) {
  x <- x[is.finite(x)]
  n <- length(x)
  est <- stats::sd(x)
  alpha <- (1 - conf) / 2
  df <- n - 1
  tibble(
    n = n,
    MAR = stats::median(abs(x), na.rm = TRUE),
    SDNR = est,
    LCI = sqrt(df * est^2 / stats::qchisq(1 - alpha, df)),
    HCI = sqrt(df * est^2 / stats::qchisq(alpha, df))
  )
}

repo_dir <- Sys.getenv("OPAL_DOCUMENTATION_DIR", "/home/darcy/Projects/opal-documentation")
opal_dir <- Sys.getenv("OPAL_SOURCE_DIR", "/home/darcy/Projects/opal")
cache_dir <- file.path(repo_dir, "assets", "cached", "bet")
out_dir <- file.path(cache_dir, "bet_cpue_tau_zero")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

base_fit <- readRDS(file.path(cache_dir, "opal_bet_fit.rds"))
if (!nzchar(Sys.getenv("OPAL_BET_LF_LOG_B0_LOWER", unset = ""))) {
  Sys.setenv(OPAL_BET_LF_LOG_B0_LOWER = "12")
}

setwd(opal_dir)
msg("Sourcing opal BET diagnostic helpers")
sys.source("inst/scripts/bet_wf_index_diagnostics.R", globalenv())
combo$helper$source_opal_r()
inputs <- combo$helper$load_bet_inputs()
load("data/wcpo_bet_wf.rda")

lf_fisheries <- c(8L, 9L, 10L, 11L, 12L, 13L, 14L)
wf_fisheries <- c(1L, 2L, 3L, 4L, 6L, 7L, 15L)
lf_var_adjust_scalar <- 400
lf_var_adjust_override <- c(`8` = 25000, `10` = 1600)
wf_var_adjust_scalar <- 10000
wf_var_adjust_override <- c(`3` = 40000, `6` = 40000, `7` = 40000)
tau_value <- as.numeric(Sys.getenv("OPAL_BET_CPUE_TAU", "1e-8"))

msg("Preparing retained BET data with CPUE tau=", tau_value)
data <- combo$helper$add_length_bins(inputs$data)
data <- combo$helper$set_quarter_cpue(data)

lf_wide <- combo$helper$make_lf_wide(inputs$lf)
lf_var_adjust <- rep(lf_var_adjust_scalar, data$n_fishery)
lf_var_adjust[as.integer(names(lf_var_adjust_override))] <- as.numeric(lf_var_adjust_override)
data <- prep_lf_data(
  data = data,
  lf_wide = lf_wide,
  lf_keep_fisheries = lf_fisheries,
  lf_var_adjust = lf_var_adjust,
  lf_switch = 1L
)

wf_wide <- make_wf_wide(wcpo_bet_wf)
data$wt_bin_start <- 1
data$wt_bin_width <- 1
data$n_wt <- 200L
wf_var_adjust <- rep(wf_var_adjust_scalar, data$n_fishery)
wf_var_adjust[as.integer(names(wf_var_adjust_override))] <- as.numeric(wf_var_adjust_override)
data <- prep_wf_data(
  data = data,
  wf_wide = wf_wide,
  wf_keep_fisheries = wf_fisheries,
  wf_switch = 1L,
  wf_var_adjust = wf_var_adjust
)
data$lf_switch <- 1L
data$wf_switch <- 1L
data$priors <- list()

parameters <- combo$helper$make_parameters(data, inputs$parameters)
parameters$log_cpue_tau <- rep(log(tau_value), data$n_index)

map_info <- make_index_map(
  data = data,
  parameters = parameters,
  estimate_index_sel = TRUE,
  index_sel_cols = c(1L),
  sel_specs = base_fit$diagnostics$selectivity_raw$sel_specs
)

obj <- MakeADFun(
  func = cmb(opal_model, data),
  parameters = parameters,
  map = map_info$map
)
obj$env$tracemgc <- FALSE
bounds <- combo$helper$make_bounds(obj, data, map_info$map, map_info$map_sel)

start <- obj$par
if (length(base_fit$opt$par) == length(start) &&
    identical(names(base_fit$opt$par), names(start))) {
  start <- base_fit$opt$par
  msg("Warm-starting from retained BET cache")
}

msg("Start nll=", signif(obj$fn(start), 8), ", npar=", length(start))
control <- list(
  eval.max = as.integer(Sys.getenv("OPAL_BET_CPUE_TAU_EVAL_MAX", "2000")),
  iter.max = as.integer(Sys.getenv("OPAL_BET_CPUE_TAU_ITER_MAX", "2000"))
)
max_restarts <- as.integer(Sys.getenv("OPAL_BET_CPUE_TAU_MAX_RESTARTS", "4"))

runtime <- system.time({
  opt <- list(par = start, convergence = NA_integer_, message = "not run")
  for (i in seq_len(max_restarts)) {
    opt <- nlminb(
      start = opt$par,
      objective = obj$fn,
      gradient = obj$gr,
      lower = bounds$lower,
      upper = bounds$upper,
      control = control
    )
    grad <- obj$gr(opt$par)
    raw_grad <- max(abs(grad), na.rm = TRUE)
    kkt_grad <- max_kkt_gradient(opt$par, grad, bounds$lower, bounds$upper)
    msg(
      "restart ", i,
      ": convergence=", opt$convergence,
      ", nll=", signif(opt$objective, 8),
      ", raw_grad=", signif(raw_grad, 6),
      ", kkt_grad=", signif(kkt_grad, 6)
    )
    if (is.finite(kkt_grad) && kkt_grad < 1e-2) break
  }
})

obj$env$last.par.best <- opt$par
obj$env$last.par <- opt$par
report <- obj$report(opt$par)
grad <- obj$gr(opt$par)
kkt_grad <- max_kkt_gradient(opt$par, grad, bounds$lower, bounds$upper)
cpue_residual <- (log(data$cpue_data$value) - log(as.numeric(report$cpue_pred))) /
  as.numeric(report$cpue_sigma)
cpue_sdnr <- bind_rows(
  sdnr_summary(cpue_residual) |> mutate(index = "All indices", .before = 1),
  bind_rows(lapply(sort(unique(data$cpue_data$month)), function(m) {
    rows <- data$cpue_data$month == m
    sdnr_summary(cpue_residual[rows]) |>
      mutate(index = month.name[m], .before = 1)
  }))
)

summary <- tibble(
  label = "cpue_tau_zero",
  tau = tau_value,
  convergence = opt$convergence,
  message = opt$message,
  objective = opt$objective,
  raw_max_gradient = max(abs(grad), na.rm = TRUE),
  kkt_max_gradient = kkt_grad,
  n_par = length(opt$par),
  runtime_seconds = unname(runtime[["elapsed"]]),
  B0 = as.numeric(report$B0),
  final_spawning_biomass = tail(as.numeric(report$spawning_biomass_y), 1),
  final_static_depletion = tail(as.numeric(report$static_depletion_y), 1),
  final_dynamic_depletion = tail(as.numeric(report$dynamic_depletion_y), 1),
  lp_prior = component_sum(report, "lp_prior"),
  lp_penalty = component_sum(report, "lp_penalty"),
  lp_rec = component_sum(report, "lp_rec"),
  lp_cpue = component_sum(report, "lp_cpue"),
  lp_lf = component_sum(report, "lp_lf"),
  lp_wf = component_sum(report, "lp_wf"),
  cpue_log_rmse = sqrt(mean((log(data$cpue_data$value) -
                               log(as.numeric(report$cpue_pred)))^2,
                             na.rm = TRUE)),
  cpue_sdnr_all = cpue_sdnr$SDNR[cpue_sdnr$index == "All indices"],
  cpue_sigma_mean = mean(as.numeric(report$cpue_sigma), na.rm = TRUE)
)

result <- list(
  label = "cpue_tau_zero",
  tau = tau_value,
  data = data,
  parameters = parameters,
  map = map_info$map,
  map_sel = map_info$map_sel,
  obj = obj,
  bounds = bounds,
  opt = opt,
  report = report,
  raw_grad = max(abs(grad), na.rm = TRUE),
  kkt_grad = kkt_grad,
  runtime = runtime,
  cpue_sdnr = cpue_sdnr,
  summary = summary
)

saveRDS(result, file.path(out_dir, "tau_zero_result.rds"))
write.csv(summary, file.path(out_dir, "tau_zero_summary.csv"), row.names = FALSE)
write.csv(cpue_sdnr, file.path(out_dir, "tau_zero_cpue_sdnr.csv"), row.names = FALSE)

msg("Wrote ", out_dir)
print(summary)
print(cpue_sdnr)
