#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(RTMB)
  library(dplyr)
  library(tidyr)
  library(purrr)
})

opal_dir <- Sys.getenv("OPAL_SOURCE_DIR", unset = "/home/darcy/Projects/opal")
out_dir <- file.path("assets", "cached", "bet")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

for (file in list.files(file.path(opal_dir, "R"), pattern = "\\.R$", full.names = TRUE)) {
  sys.source(file, envir = .GlobalEnv)
}

load(file.path(opal_dir, "data", "wcpo_bet_data.rda"))
load(file.path(opal_dir, "data", "wcpo_bet_parameters.rda"))
load(file.path(opal_dir, "data", "wcpo_bet_lf.rda"))
load(file.path(opal_dir, "data", "wcpo_bet_wf.rda"))

fit <- readRDS(file.path(out_dir, "opal_bet_fit.rds"))
cpue_tau <- as.numeric(Sys.getenv("OPAL_BET_CPUE_TAU", "1e-8"))
if (!is.finite(cpue_tau) || cpue_tau <= 0) {
  stop("OPAL_BET_CPUE_TAU must be a positive finite value")
}
xmodel_biomass <- readRDS(file.path(out_dir, "wcpo_bet_xmodel_biomass.rds"))

add_length_bins <- function(data) {
  data$len_lower <- seq(data$len_bin_start, by = data$len_bin_width, length.out = data$n_len)
  data$len_upper <- data$len_lower + data$len_bin_width
  data$len_mid <- data$len_lower + data$len_bin_width / 2
  data
}

set_quarter_cpue <- function(data) {
  data$cpue_data$index <- as.integer(factor(data$cpue_data$month, levels = c(2, 5, 8, 11)))
  data$n_index <- 4L
  data
}

make_lf_wide <- function(lf) {
  lf |>
    pivot_wider(
      id_cols = c(fishery, year, month, ts),
      names_from = bin,
      values_from = value,
      values_fill = 0
    ) |>
    arrange(fishery, ts)
}

make_wf_wide <- function(wf) {
  wf |>
    pivot_wider(
      id_cols = c(fishery, year, month, ts),
      names_from = bin,
      values_from = value,
      values_fill = 0
    ) |>
    arrange(fishery, ts)
}

prepare_full_data <- function(lf_switch = 1L, wf_switch = 1L, sex_ratio = NULL) {
  data <- wcpo_bet_data |>
    add_length_bins() |>
    set_quarter_cpue()

  lf_var_adjust <- rep(80, data$n_fishery)
  lf_var_adjust[8] <- 5000
  data <- prep_lf_data(
    data = data,
    lf_wide = make_lf_wide(wcpo_bet_lf),
    lf_keep_fisheries = c(8, 9, 10, 11, 12, 13, 14),
    lf_var_adjust = lf_var_adjust,
    lf_switch = 1L
  )

  data$wt_bin_start <- 1
  data$wt_bin_width <- 1
  data$n_wt <- 200L
  data <- prep_wf_data(
    data = data,
    wf_wide = make_wf_wide(wcpo_bet_wf),
    wf_keep_fisheries = c(1, 2, 3, 4, 6, 7, 15),
    wf_switch = 1L,
    wf_var_adjust = rep(2000, data$n_fishery)
  )

  data$lf_switch <- as.integer(lf_switch)
  data$wf_switch <- as.integer(wf_switch)
  data$priors <- list()
  if (!is.null(sex_ratio)) data$sex_ratio <- rep(sex_ratio, data$n_age)
  data
}

make_fitted_parameters <- function(data, par_sel = fit$diagnostics$selectivity_raw$par_est) {
  p <- fit$diagnostics$selectivity_raw$opt_par
  list(
    log_B0 = unname(p[names(p) == "log_B0"]),
    log_h = as.numeric(wcpo_bet_parameters$log_h),
    log_sigma_r = as.numeric(wcpo_bet_parameters$log_sigma_r),
    log_cpue_q = unname(p[names(p) == "log_cpue_q"]),
    cpue_creep = rep(as.numeric(wcpo_bet_parameters$cpue_creep), data$n_index),
    log_cpue_tau = rep(log(cpue_tau), data$n_index),
    log_cpue_omega = rep(as.numeric(wcpo_bet_parameters$log_cpue_omega), data$n_index),
    log_lf_tau = as.numeric(log(rep(0.1, data$n_fishery))),
    log_wf_tau = rep(0, data$n_fishery),
    log_L1 = as.numeric(wcpo_bet_parameters$log_L1),
    log_L2 = as.numeric(wcpo_bet_parameters$log_L2),
    log_k = as.numeric(wcpo_bet_parameters$log_k),
    log_CV1 = as.numeric(wcpo_bet_parameters$log_CV1),
    log_CV2 = as.numeric(wcpo_bet_parameters$log_CV2),
    par_sel = par_sel,
    rdev_y = unname(p[names(p) == "rdev_y"])
  )
}

make_fitted_map <- function(data, par_sel_map = fit$diagnostics$selectivity_raw$map_sel) {
  list(
    log_h = factor(NA),
    log_sigma_r = factor(NA),
    log_cpue_q = factor(seq_len(data$n_index)),
    cpue_creep = factor(rep(NA, data$n_index)),
    log_cpue_tau = factor(rep(NA, data$n_index)),
    log_cpue_omega = factor(rep(NA, data$n_index)),
    log_lf_tau = factor(rep(NA_integer_, data$n_fishery)),
    log_wf_tau = factor(rep(NA_integer_, data$n_fishery)),
    log_L1 = factor(NA),
    log_L2 = factor(NA),
    log_k = factor(NA),
    log_CV1 = factor(NA),
    log_CV2 = factor(NA),
    par_sel = factor(par_sel_map),
    rdev_y = factor(seq_along(fit$diagnostics$selectivity_raw$opt_par[names(fit$diagnostics$selectivity_raw$opt_par) == "rdev_y"]))
  )
}

component_row <- function(label, obj, par, data, optimize_q = FALSE) {
  if (optimize_q) {
    q_pos <- which(names(par) == "log_cpue_q")
    fixed_par <- par
    q_opt <- nlminb(
      start = par[q_pos],
      objective = function(q) {
        p <- fixed_par
        p[q_pos] <- q
        obj$fn(p)
      },
      lower = rep(log(0.1), length(q_pos)),
      upper = rep(log(10), length(q_pos)),
      control = list(eval.max = 100, iter.max = 100)
    )
    par[q_pos] <- q_opt$par
  } else {
    q_opt <- NULL
  }

  rep <- obj$report(par)
  tibble(
    label = label,
    optimize_q = optimize_q,
    objective = obj$fn(par),
    B0 = as.numeric(rep$B0),
    log_B0 = log(as.numeric(rep$B0)),
    final_spawning_biomass = tail(as.numeric(rep$spawning_biomass_y), 1),
    final_static_depletion = tail(as.numeric(rep$static_depletion_y), 1),
    final_dynamic_depletion = tail(as.numeric(rep$dynamic_depletion_y), 1),
    lp_prior = sum(as.numeric(rep$lp_prior)),
    lp_penalty = sum(as.numeric(rep$lp_penalty)),
    lp_rec = sum(as.numeric(rep$lp_rec)),
    lp_cpue = sum(as.numeric(rep$lp_cpue)),
    lp_lf = sum(as.numeric(rep$lp_lf)),
    lp_wf = sum(as.numeric(rep$lp_wf)),
    cpue_log_rmse = sqrt(mean((log(data$cpue_data$value) - log(as.numeric(rep$cpue_pred)))^2)),
    max_harvest_rate = max(as.numeric(rep$hrate_ysa), na.rm = TRUE),
    q_min = min(exp(par[names(par) == "log_cpue_q"])),
    q_max = max(exp(par[names(par) == "log_cpue_q"])),
    q_convergence = if (is.null(q_opt)) NA_integer_ else q_opt$convergence
  )
}

fishery_harvest_summary <- function(rep, data) {
  out <- vector("list", data$n_fishery)
  for (f in seq_len(data$n_fishery)) {
    den <- numeric(data$n_year)
    obs <- numeric(data$n_year)
    for (y in seq_len(data$n_year)) {
      n <- rep$number_ysa[y, 1, ]
      sel <- rep$sel_fya[f, y, ]
      if (data$catch_units_f[f] == 1) {
        den[y] <- sum(n * sel * rep$weight_fya_mod[f, y, ])
      } else {
        den[y] <- sum(n * sel)
      }
      obs[y] <- sum(data$catch_obs_ysf[y, , f])
    }
    h <- ifelse(den > 0, obs / den, NA_real_)
    out[[f]] <- tibble(
      fishery = f,
      catch_units = if_else(data$catch_units_f[f] == 1, "weight", "numbers"),
      total_catch = sum(obs),
      mean_harvest_fraction = mean(h[obs > 0], na.rm = TRUE),
      max_harvest_fraction = max(h[obs > 0], na.rm = TRUE),
      active_years = sum(obs > 0)
    )
  }
  bind_rows(out)
}

selectivity_contrast <- function(data) {
  sel_at_length <- function(f, par) {
    if (fit$diagnostics$selectivity_raw$sel_type_f[f] == 1L) {
      as.numeric(sel_logistic(data$len_mid, par))
    } else {
      as.numeric(sel_double_normal(data$len_mid, par))
    }
  }
  sel_est <- map(seq_len(data$n_fishery), ~ sel_at_length(.x, fit$diagnostics$selectivity_raw$par_est[.x, ]))
  sel_ref <- map(seq_len(data$n_fishery), ~ sel_at_length(.x, fit$diagnostics$selectivity_raw$par_ref[.x, ]))
  mature_len <- data$maturity > 0.5
  tibble(
    fishery = seq_len(data$n_fishery),
    mean_sel_est = map_dbl(sel_est, mean),
    mean_sel_ref = map_dbl(sel_ref, mean),
    mature_len_sel_est = map_dbl(sel_est, ~ mean(.x[mature_len])),
    mature_len_sel_ref = map_dbl(sel_ref, ~ mean(.x[mature_len])),
    high_len_sel_est = map_dbl(sel_est, ~ mean(.x[data$len_mid >= 120])),
    high_len_sel_ref = map_dbl(sel_ref, ~ mean(.x[data$len_mid >= 120]))
  )
}

cat("Preparing retained BET data and AD object...\n")
data <- prepare_full_data()
pars <- make_fitted_parameters(data)
map <- make_fitted_map(data)
obj <- MakeADFun(func = cmb(opal_model, data), parameters = pars, map = map)
base_par <- obj$par
base_rep <- obj$report(base_par)

xmodel_b0 <- xmodel_biomass |>
  mutate(B0 = ssb / depletion) |>
  group_by(model) |>
  summarise(B0 = first(B0), .groups = "drop")

targets <- bind_rows(
  tibble(label = "opal_cached", B0 = as.numeric(fit$report$B0)),
  xmodel_b0 |> transmute(label = paste0("xmodel_", model), B0)
) |>
  distinct(label, .keep_all = TRUE)

profile_rows <- map_dfr(seq_len(nrow(targets)), function(i) {
  p <- base_par
  p[names(p) == "log_B0"] <- log(targets$B0[i])
  bind_rows(
    component_row(targets$label[i], obj, p, data, optimize_q = FALSE),
    component_row(targets$label[i], obj, p, data, optimize_q = TRUE)
  )
})

if (identical(Sys.getenv("OPAL_BET_RUN_SWITCHES", unset = "false"), "true")) {
  switch_rows <- map_dfr(
    list(
      retained = c(1L, 1L),
      no_wf = c(1L, 0L),
      no_lf = c(0L, 1L),
      cpue_catch_only = c(0L, 0L)
    ),
    function(sw) {
      d <- prepare_full_data(lf_switch = sw[[1]], wf_switch = sw[[2]])
      o <- MakeADFun(func = cmb(opal_model, d), parameters = make_fitted_parameters(d), map = make_fitted_map(d))
      component_row(paste0("switch_lf", sw[[1]], "_wf", sw[[2]]), o, o$par, d, optimize_q = FALSE)
    },
    .id = "switch"
  )
} else {
  switch_rows <- tibble(
    switch = character(),
    label = character(),
    optimize_q = logical(),
    objective = double(),
    B0 = double(),
    log_B0 = double(),
    final_spawning_biomass = double(),
    final_static_depletion = double(),
    final_dynamic_depletion = double(),
    lp_prior = double(),
    lp_penalty = double(),
    lp_rec = double(),
    lp_cpue = double(),
    lp_lf = double(),
    lp_wf = double(),
    cpue_log_rmse = double(),
    max_harvest_rate = double(),
    q_min = double(),
    q_max = double(),
    q_convergence = integer()
  )
}

harvest_summary <- fishery_harvest_summary(base_rep, data)
sel_contrast <- selectivity_contrast(data)

out <- list(
  created = Sys.time(),
  retained_objective = obj$fn(base_par),
  cached_objective = fit$opt$objective,
  xmodel_b0 = xmodel_b0,
  profile = profile_rows,
  likelihood_switches = switch_rows,
  harvest_summary = harvest_summary,
  selectivity_contrast = sel_contrast
)

saveRDS(out, file.path(out_dir, "opal_bet_biomass_scale_checks.rds"))
write.csv(profile_rows, file.path(out_dir, "opal_bet_biomass_scale_profile.csv"), row.names = FALSE)
write.csv(switch_rows, file.path(out_dir, "opal_bet_biomass_scale_likelihood_switches.csv"), row.names = FALSE)
write.csv(harvest_summary, file.path(out_dir, "opal_bet_biomass_scale_harvest.csv"), row.names = FALSE)
write.csv(sel_contrast, file.path(out_dir, "opal_bet_biomass_scale_selectivity.csv"), row.names = FALSE)

cat("\nB0 profile/evaluation summary:\n")
print(profile_rows |>
  select(label, optimize_q, objective, B0, final_static_depletion,
         lp_penalty, lp_cpue, lp_lf, lp_wf, cpue_log_rmse, max_harvest_rate,
         q_min, q_max))

if (nrow(switch_rows)) {
  cat("\nLikelihood-switch summary at retained fitted parameters:\n")
  print(switch_rows |>
    select(switch, objective, B0, lp_penalty, lp_cpue, lp_lf, lp_wf, cpue_log_rmse))
}

cat("\nHarvest summary:\n")
print(harvest_summary)
