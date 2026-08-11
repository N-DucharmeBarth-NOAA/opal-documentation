#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(RTMB)
  library(dplyr)
  library(tidyr)
  library(purrr)
})

msg <- function(...) {
  cat(format(Sys.time(), "%H:%M:%S"), ..., "\n")
  flush(stdout())
}

repo_dir <- "/home/darcy/Projects/opal-documentation"
opal_dir <- Sys.getenv("OPAL_SOURCE_DIR", "/home/darcy/Projects/opal")
out_dir <- Sys.getenv(
  "OPAL_BET_DW_OUT_DIR",
  file.path(repo_dir, "assets", "cached", "bet", "bet_lfwf_downweight10_raw")
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

mult <- as.numeric(Sys.getenv("OPAL_BET_LFWF_DOWNWEIGHT_MULT", "10"))
if (!is.finite(mult) || mult <= 0) {
  stop("OPAL_BET_LFWF_DOWNWEIGHT_MULT must be positive")
}

if (!nzchar(Sys.getenv("OPAL_BET_LF_LOG_B0_LOWER", unset = ""))) {
  Sys.setenv(OPAL_BET_LF_LOG_B0_LOWER = "12")
}

cache_path <- file.path(repo_dir, "assets", "cached", "bet", "opal_bet_fit.rds")
cached_fit <- readRDS(cache_path)

setwd(opal_dir)
sys.source("inst/scripts/bet_wf_index_diagnostics.R", globalenv())
combo$helper$source_opal_r()
inputs <- combo$helper$load_bet_inputs()
load("data/wcpo_bet_wf.rda")

wf_wide <- make_wf_wide(wcpo_bet_wf)
lf_wide <- combo$helper$make_lf_wide(inputs$lf)
sel_specs <- cached_fit$diagnostics$selectivity_raw$sel_specs

base_prepare_lf_wf_data <- prepare_lf_wf_data
prepare_lf_wf_data <- function(data, lf_wide, wf_wide,
                               lf_fisheries = c(8L, 9L, 10L, 11L, 12L, 13L, 14L),
                               wf_fisheries = c(1L, 2L, 3L, 4L, 6L, 7L, 15L),
                               lf_var_adjust_scalar = 80 * mult,
                               lf_var_adjust_override = c(`8` = 5000 * mult),
                               wf_var_adjust_scalar = 2000 * mult) {
  base_prepare_lf_wf_data(
    data = data,
    lf_wide = lf_wide,
    wf_wide = wf_wide,
    lf_fisheries = lf_fisheries,
    wf_fisheries = wf_fisheries,
    lf_var_adjust_scalar = lf_var_adjust_scalar,
    lf_var_adjust_override = lf_var_adjust_override,
    wf_var_adjust_scalar = wf_var_adjust_scalar
  )
}

msg("Building 15-fishery BET fit with LF/WF downweight multiplier ", mult)
msg("LF adjust=", 80 * mult, ", LF8 adjust=", 5000 * mult,
    ", WF adjust=", 2000 * mult)

built <- build_fit(
  base_data = inputs$data,
  wcpo_pars = inputs$parameters,
  lf_wide = lf_wide,
  wf_wide = wf_wide,
  wf_var_adjust_scalar = 2000 * mult,
  lf_fisheries = c(8L, 9L, 10L, 11L, 12L, 13L, 14L),
  wf_fisheries = c(1L, 2L, 3L, 4L, 6L, 7L, 15L),
  lf_var_adjust_override = c(`8` = 5000 * mult),
  sel_specs = sel_specs
)

obj <- built$obj
start <- obj$par
if (length(cached_fit$opt$par) == length(start) &&
    identical(names(cached_fit$opt$par), names(start))) {
  start <- cached_fit$opt$par
  msg("Warm-starting from retained BET cache")
} else {
  msg("Cached parameter vector does not match; using configured starts")
}

control <- list(
  eval.max = as.integer(Sys.getenv("OPAL_BET_DW_EVAL_MAX", "3000")),
  iter.max = as.integer(Sys.getenv("OPAL_BET_DW_ITER_MAX", "3000"))
)
max_restarts <- as.integer(Sys.getenv("OPAL_BET_DW_MAX_RESTARTS", "6"))

msg("Start nll=", signif(obj$fn(start), 8), ", npar=", length(start))
runtime <- system.time({
  opt <- list(par = start, convergence = NA_integer_, message = "not run")
  for (i in seq_len(max_restarts)) {
    opt <- nlminb(
      start = opt$par,
      objective = obj$fn,
      gradient = obj$gr,
      lower = built$bounds$lower,
      upper = built$bounds$upper,
      control = control
    )
    grad <- obj$gr(opt$par)
    raw_grad <- max(abs(grad))
    kkt_grad <- combo$max_kkt_gradient(
      opt$par, grad, built$bounds$lower, built$bounds$upper
    )
    msg("restart ", i,
        ": convergence=", opt$convergence,
        ", nll=", signif(opt$objective, 8),
        ", raw_grad=", signif(raw_grad, 6),
        ", kkt_grad=", signif(kkt_grad, 6))
    if (is.finite(kkt_grad) && kkt_grad < 1e-2) break
  }
})

obj$env$last.par.best <- opt$par
obj$env$last.par <- opt$par
rep <- obj$report(opt$par)
par_list <- obj$env$parList(opt$par)
grad <- obj$gr(opt$par)
kkt_grad <- combo$max_kkt_gradient(
  opt$par, grad, built$bounds$lower, built$bounds$upper
)

result <- list(
  label = paste0("lfwf_downweight", format(mult, trim = TRUE)),
  lf_var_adjust_scalar = 80 * mult,
  wf_var_adjust = 2000 * mult,
  lf_var_adjust_override = c(`8` = 5000 * mult),
  estimate_index_sel = TRUE,
  index_sel_cols = c(1L),
  lf_fisheries = c(8L, 9L, 10L, 11L, 12L, 13L, 14L),
  wf_fisheries = c(1L, 2L, 3L, 4L, 6L, 7L, 15L),
  sel_specs = sel_specs,
  fit = built,
  opt = opt,
  report = rep,
  par_list = par_list,
  raw_grad = max(abs(grad)),
  kkt_grad = kkt_grad,
  runtime = runtime,
  downweight_multiplier = mult
)

saveRDS(list(result), file.path(out_dir, "wf_index_results.rds"))
metrics <- make_metrics(result) %>%
  mutate(
    lf_var_adjust_scalar = 80 * mult,
    lf8_var_adjust = 5000 * mult,
    runtime_seconds = unname(runtime[["elapsed"]])
  )
write.csv(metrics, file.path(out_dir, "wf_index_metrics.csv"), row.names = FALSE)
write.csv(make_lf_metrics(result), file.path(out_dir, "wf_index_lf_metrics.csv"), row.names = FALSE)
write.csv(make_wf_metrics(result), file.path(out_dir, "wf_index_wf_metrics.csv"), row.names = FALSE)

msg("Wrote ", file.path(out_dir, "wf_index_results.rds"))
print(metrics %>%
  select(label, n_par, n_lf_obs, n_wf_obs, convergence, nll,
         raw_max_gradient, kkt_max_gradient, B0, lp_cpue, lp_lf, lp_wf,
         log_cpue_rmse, final_spawning_biomass_ratio, runtime_seconds))
