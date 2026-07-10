#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(RTMB)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
})

msg <- function(...) {
  cat(format(Sys.time(), "%H:%M:%S"), ..., "\n")
  flush(stdout())
}

out_dir <- Sys.getenv("OPAL_BET_AGG6_OUT_DIR", "/tmp/opal_bet_agg6")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

setwd("/home/darcy/Projects/opal")
sys.source("inst/scripts/bet_wf_index_diagnostics.R", globalenv())
combo$helper$source_opal_r()
if (!nzchar(Sys.getenv("OPAL_BET_LF_LOG_B0_LOWER", unset = ""))) {
  Sys.setenv(OPAL_BET_LF_LOG_B0_LOWER = "12")
}

inputs <- combo$helper$load_bet_inputs()
load("data/wcpo_bet_wf.rda")
cpue_tau <- as.numeric(Sys.getenv("OPAL_BET_CPUE_TAU", "1e-8"))
if (!is.finite(cpue_tau) || cpue_tau <= 0) {
  stop("OPAL_BET_CPUE_TAU must be a positive finite value")
}

base_old <- inputs$data
pars_old <- inputs$parameters
lf_old <- inputs$lf
wf_old <- wcpo_bet_wf

groups <- list(
  `1` = 8L,
  `2` = 10L,
  `3` = c(9L, 12L, 13L, 14L),
  `4` = 11L,
  `5` = 1:7,
  `6` = 15L
)

old_to_new <- rep(NA_integer_, base_old$n_fishery)
for (g in names(groups)) old_to_new[groups[[g]]] <- as.integer(g)
stopifnot(!anyNA(old_to_new))

catch_old_f <- as.numeric(apply(base_old$catch_obs_ysf, 3, sum))
group_summary <- imap_dfr(groups, function(old_f, g_chr) {
  g <- as.integer(g_chr)
  tibble(
    group = g,
    old_fisheries = paste(old_f, collapse = ","),
    catch_units = unique(base_old$catch_units_f[old_f]),
    sel_type = if_else(unique(base_old$sel_type_f[old_f]) == 1L, "logistic", "double-normal"),
    total_catch = sum(catch_old_f[old_f])
  )
}) |>
  group_by(catch_units) |>
  mutate(percent = 100 * total_catch / sum(total_catch)) |>
  ungroup()

if (any(lengths(strsplit(group_summary$sel_type, ",", fixed = TRUE)) > 1L)) {
  stop("A group has mixed selectivity types")
}
if (any(lengths(strsplit(as.character(group_summary$catch_units), ",", fixed = TRUE)) > 1L)) {
  stop("A group has mixed catch units")
}

msg("six-fishery grouping")
print(group_summary)
write.csv(group_summary, file.path(out_dir, "grouping_summary.csv"), row.names = FALSE)

aggregate_catch <- function(catch_ysf, old_to_new, n_new) {
  dims <- dim(catch_ysf)
  out <- array(0, dim = c(dims[1], dims[2], n_new))
  for (old_f in seq_along(old_to_new)) {
    new_f <- old_to_new[[old_f]]
    out[, , new_f] <- out[, , new_f] + catch_ysf[, , old_f]
  }
  out
}

aggregate_comp <- function(x, old_to_new) {
  x |>
    mutate(
      old_fishery = fishery,
      fishery = old_to_new[old_fishery]
    ) |>
    group_by(fishery, year, month, ts, bin) |>
    summarise(value = sum(value), .groups = "drop") |>
    arrange(fishery, ts, bin)
}

aggregate_par_sel <- function(par_sel, groups, weights) {
  out <- matrix(NA_real_, nrow = length(groups), ncol = ncol(par_sel))
  for (g_chr in names(groups)) {
    g <- as.integer(g_chr)
    old_f <- groups[[g_chr]]
    w <- weights[old_f]
    if (!sum(w) > 0) w <- rep(1, length(old_f))
    w <- w / sum(w)
    out[g, ] <- colSums(par_sel[old_f, , drop = FALSE] * w)
  }
  colnames(out) <- colnames(par_sel)
  out
}

base_new <- base_old
base_new$n_fishery <- length(groups)
base_new$catch_obs_ysf <- aggregate_catch(base_old$catch_obs_ysf, old_to_new, length(groups))
base_new$catch_units_f <- as.numeric(group_summary$catch_units)
base_new$sel_type_f <- c(2L, 2L, 2L, 1L, 2L, 1L)
base_new$cpue_data <- base_new$cpue_data |>
  mutate(fishery = old_to_new[fishery])

pars_new <- pars_old
pars_new$par_sel <- aggregate_par_sel(pars_old$par_sel, groups, catch_old_f)

lf_new <- aggregate_comp(lf_old, old_to_new)
wf_new <- aggregate_comp(wf_old, old_to_new)
lf_wide_new <- combo$helper$make_lf_wide(lf_new)
wf_wide_new <- make_wf_wide(wf_new)

comp_upweight <- as.numeric(Sys.getenv("OPAL_BET_AGG6_COMP_UPWEIGHT", "1"))
if (!is.finite(comp_upweight) || comp_upweight <= 0) {
  stop("OPAL_BET_AGG6_COMP_UPWEIGHT must be a positive number")
}
agg6_lf_var_adjust_scalar <- as.numeric(Sys.getenv(
  "OPAL_BET_AGG6_LF_VAR_ADJUST",
  as.character(80 / comp_upweight)
))
agg6_lf_group1_var_adjust <- as.numeric(Sys.getenv(
  "OPAL_BET_AGG6_LF_GROUP1_VAR_ADJUST",
  as.character(5000 / comp_upweight)
))
agg6_wf_var_adjust_scalar <- as.numeric(Sys.getenv(
  "OPAL_BET_AGG6_WF_VAR_ADJUST",
  as.character(2000 / comp_upweight)
))

base_prepare_lf_wf_data <- prepare_lf_wf_data
prepare_lf_wf_data <- function(data, lf_wide, wf_wide,
                               lf_fisheries = c(9L, 10L, 14L),
                               wf_fisheries = 15L,
                               lf_var_adjust_scalar = agg6_lf_var_adjust_scalar,
                               lf_var_adjust_override = numeric(),
                               wf_var_adjust_scalar = agg6_wf_var_adjust_scalar) {
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

variant <- Sys.getenv("OPAL_BET_AGG6_VARIANT", "full_g5")
sel_specs <- switch(
  variant,
  g5_no_desc = list(
    `1` = 1L,
    `2` = c(1L, 4L),
    `3` = c(1L, 4L),
    `4` = c(1L, 2L),
    `5` = c(1L, 3L),
    `6` = 1L
  ),
  g5_peak_only = list(
    `1` = 1L,
    `2` = c(1L, 4L),
    `3` = c(1L, 4L),
    `4` = c(1L, 2L),
    `5` = 1L,
    `6` = 1L
  ),
  list(
    `1` = 1L,
    `2` = c(1L, 4L),
    `3` = c(1L, 4L),
    `4` = c(1L, 2L),
    `5` = c(1L, 3L, 4L),
    `6` = 1L
  )
)
fit_label <- Sys.getenv(
  "OPAL_BET_AGG6_LABEL",
  paste0("agg6_", variant, "_upw", format(comp_upweight, trim = TRUE), "_wf", agg6_wf_var_adjust_scalar)
)

msg("aggregated composition rows: LF=", nrow(lf_wide_new), ", WF=", nrow(wf_wide_new))
msg("composition weighting: upweight=", comp_upweight,
    ", lf_adjust=", agg6_lf_var_adjust_scalar,
    ", lf_group1_adjust=", agg6_lf_group1_var_adjust,
    ", wf_adjust=", agg6_wf_var_adjust_scalar,
    ", cpue_tau=", cpue_tau)
msg("fitting aggregated BET model")

fit_args <- list(
  base_data = base_new,
  wcpo_pars = pars_new,
  lf_wide = lf_wide_new,
  wf_wide = wf_wide_new,
  wf_var_adjust_scalar = agg6_wf_var_adjust_scalar,
  lf_fisheries = 1:4,
  wf_fisheries = 5:6,
  lf_var_adjust_override = c(`1` = agg6_lf_group1_var_adjust),
  sel_specs = sel_specs
)

build_fit_switch <- function(fit_args, cpue_switch = 1L, lf_switch = 1L, wf_switch = 1L) {
  data <- prepare_lf_wf_data(
    data = fit_args$base_data,
    lf_wide = fit_args$lf_wide,
    wf_wide = fit_args$wf_wide,
    lf_fisheries = fit_args$lf_fisheries,
    wf_fisheries = fit_args$wf_fisheries,
    lf_var_adjust_override = fit_args$lf_var_adjust_override,
    wf_var_adjust_scalar = fit_args$wf_var_adjust_scalar
  )
  data$cpue_switch <- as.integer(cpue_switch)
  data$lf_switch <- as.integer(lf_switch)
  data$wf_switch <- as.integer(wf_switch)
  parameters <- combo$helper$make_parameters(data, fit_args$wcpo_pars)
  parameters$log_cpue_tau <- rep(log(cpue_tau), data$n_index)
  map_info <- make_index_map(
    data = data,
    parameters = parameters,
    sel_specs = fit_args$sel_specs
  )
  obj <- MakeADFun(
    func = cmb(opal_model, data),
    parameters = parameters,
    map = map_info$map
  )
  obj$env$tracemgc <- FALSE
  bounds <- combo$helper$make_bounds(obj, data, map_info$map, map_info$map_sel)
  list(data = data, parameters = parameters, map = map_info$map,
       map_sel = map_info$map_sel, obj = obj, bounds = bounds)
}

preflight <- build_fit_switch(fit_args)
start_nll <- preflight$obj$fn(preflight$obj$par)
msg("preflight start nll=", signif(start_nll, 8))
if (!is.finite(start_nll) || identical(Sys.getenv("OPAL_BET_AGG6_DEBUG"), "true")) {
  switch_checks <- tribble(
    ~label, ~cpue_switch, ~lf_switch, ~wf_switch,
    "penalty_rec_only", 0L, 0L, 0L,
    "cpue_only", 1L, 0L, 0L,
    "lf_only", 0L, 1L, 0L,
    "wf_only", 0L, 0L, 1L,
    "lf_wf_only", 0L, 1L, 1L
  ) |>
    mutate(
      nll = pmap_dbl(
        list(cpue_switch, lf_switch, wf_switch),
        function(cpue_switch, lf_switch, wf_switch) {
          b <- build_fit_switch(
            fit_args,
            cpue_switch = cpue_switch,
            lf_switch = lf_switch,
            wf_switch = wf_switch
          )
          tryCatch(b$obj$fn(b$obj$par), error = function(e) NaN)
        }
      )
    )
  print(switch_checks)
  rep0 <- tryCatch(preflight$obj$report(preflight$obj$par), error = identity)
  if (inherits(rep0, "condition")) {
    msg("preflight report failed: ", conditionMessage(rep0))
  } else {
    component_summary <- tibble(
      component = c("lp_prior", "lp_penalty", "lp_rec", "lp_cpue", "lp_lf", "lp_wf"),
      value = c(
        sum(as.numeric(rep0$lp_prior)),
        sum(as.numeric(rep0$lp_penalty)),
        sum(as.numeric(rep0$lp_rec)),
        sum(as.numeric(rep0$lp_cpue)),
        sum(as.numeric(rep0$lp_lf)),
        sum(as.numeric(rep0$lp_wf))
      ),
      n_nonfinite = c(
        sum(!is.finite(as.numeric(rep0$lp_prior))),
        sum(!is.finite(as.numeric(rep0$lp_penalty))),
        sum(!is.finite(as.numeric(rep0$lp_rec))),
        sum(!is.finite(as.numeric(rep0$lp_cpue))),
        sum(!is.finite(as.numeric(rep0$lp_lf))),
        sum(!is.finite(as.numeric(rep0$lp_wf)))
      )
    )
    print(component_summary)
    if (!is.null(rep0$catch_pred_ysf)) {
      catch_totals <- tibble(
        fishery = seq_len(dim(rep0$catch_pred_ysf)[3]),
        obs = as.numeric(apply(preflight$data$catch_obs_ysf, 3, sum)),
        pred = as.numeric(apply(rep0$catch_pred_ysf, 3, sum))
      )
      print(catch_totals)
    }
  }
  if (!is.finite(start_nll)) stop("Invalid starting objective; see preflight diagnostics above.")
}

runtime <- system.time({
  built <- preflight
  obj <- built$obj
  start <- obj$par
  start_path <- Sys.getenv("OPAL_BET_AGG6_START_RDS", "/tmp/opal_bet_agg6/agg6_result.rds")
  if (file.exists(start_path)) {
    start_result <- readRDS(start_path)$result
    if (length(start_result$opt$par) == length(start) &&
        identical(names(start_result$opt$par), names(start))) {
      start <- start_result$opt$par
      msg("warm-starting from ", start_path)
    } else {
      msg("not warm-starting; parameter vector differs from ", start_path)
    }
  }
  opt <- list(par = start, convergence = NA_integer_, message = "not run")
  control <- list(
    eval.max = as.integer(Sys.getenv("OPAL_BET_AGG6_EVAL_MAX", "1000")),
    iter.max = as.integer(Sys.getenv("OPAL_BET_AGG6_ITER_MAX", "1000"))
  )
  max_restarts <- as.integer(Sys.getenv("OPAL_BET_AGG6_MAX_RESTARTS", "4"))
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
  obj$env$last.par.best <- opt$par
  rep <- obj$report(opt$par)
  par_list <- obj$env$parList(opt$par)
  grad <- obj$gr(opt$par)
  result <- list(
    label = fit_label,
    wf_var_adjust = fit_args$wf_var_adjust_scalar,
    lf_var_adjust_override = fit_args$lf_var_adjust_override,
    estimate_index_sel = TRUE,
    index_sel_cols = c(1L, 2L),
    lf_fisheries = fit_args$lf_fisheries,
    wf_fisheries = fit_args$wf_fisheries,
    sel_specs = fit_args$sel_specs,
    fit = built,
    opt = opt,
    report = rep,
    par_list = par_list,
    raw_grad = max(abs(grad)),
    kkt_grad = combo$max_kkt_gradient(
      opt$par, grad, built$bounds$lower, built$bounds$upper
    ),
    cpue_tau = cpue_tau
  )
})

metrics <- make_metrics(result) |>
  mutate(runtime_seconds = unname(runtime[["elapsed"]]))
lf_metrics <- make_lf_metrics(result)
wf_metrics <- make_wf_metrics(result)

saveRDS(
  list(
    grouping = group_summary,
    old_to_new = old_to_new,
    result = result,
    metrics = metrics,
    lf_metrics = lf_metrics,
    wf_metrics = wf_metrics,
    runtime = runtime
  ),
  file.path(out_dir, "agg6_result.rds")
)
write.csv(metrics, file.path(out_dir, "agg6_metrics.csv"), row.names = FALSE)
write.csv(lf_metrics, file.path(out_dir, "agg6_lf_metrics.csv"), row.names = FALSE)
write.csv(wf_metrics, file.path(out_dir, "agg6_wf_metrics.csv"), row.names = FALSE)

msg("fit complete")
print(metrics |>
  select(
    label, n_par, n_lf_obs, n_wf_obs, convergence, nll,
    raw_max_gradient, kkt_max_gradient, B0,
    lp_cpue, lp_lf, lp_wf, log_cpue_rmse,
    final_spawning_biomass_ratio, runtime_seconds
  ))

print(lf_metrics |>
  select(combo, fishery, lf_prop_rmse, lf_mean_length_rmse,
         lf_mean_length_bias, bound_status, sel_est_nat))

print(wf_metrics |>
  select(label, fishery, wf_prop_rmse, wf_mean_weight_rmse,
         wf_mean_weight_bias))
