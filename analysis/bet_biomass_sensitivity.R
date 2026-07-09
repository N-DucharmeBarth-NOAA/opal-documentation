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

source_opal_r <- function(path) {
  files <- list.files(file.path(path, "R"), pattern = "\\.R$", full.names = TRUE)
  for (file in files) sys.source(file, envir = .GlobalEnv)
  invisible(files)
}

source_opal_r(opal_dir)

load(file.path(opal_dir, "data", "wcpo_bet_data.rda"))
load(file.path(opal_dir, "data", "wcpo_bet_parameters.rda"))
load(file.path(opal_dir, "data", "wcpo_bet_lf.rda"))
load(file.path(opal_dir, "data", "wcpo_bet_wf.rda"))

cached_fit <- readRDS(file.path(out_dir, "opal_bet_fit.rds"))
cpue_tau <- as.numeric(Sys.getenv("OPAL_BET_CPUE_TAU", "1e-8"))
if (!is.finite(cpue_tau) || cpue_tau <= 0) {
  stop("OPAL_BET_CPUE_TAU must be a positive finite value")
}
xmodel_biomass <- readRDS(file.path(out_dir, "wcpo_bet_xmodel_biomass.rds"))

`%||%` <- function(x, y) if (is.null(x)) y else x

add_length_bins <- function(data) {
  data$len_lower <- seq(
    from = data$len_bin_start,
    by = data$len_bin_width,
    length.out = data$n_len
  )
  data$len_upper <- data$len_lower + data$len_bin_width
  data$len_mid <- data$len_lower + data$len_bin_width / 2
  data
}

set_quarter_cpue <- function(data) {
  data$cpue_data$index <- as.integer(
    factor(data$cpue_data$month, levels = c(2, 5, 8, 11))
  )
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

cat("Preparing wide composition tables...\n")
flush(stdout())
lf_wide_global <- make_lf_wide(wcpo_bet_lf)
wf_wide_global <- make_wf_wide(wcpo_bet_wf)
cat("Wide composition tables prepared: LF rows=", nrow(lf_wide_global),
    ", WF rows=", nrow(wf_wide_global), "\n", sep = "")
flush(stdout())

make_par_sel_start <- function(data, wcpo_pars) {
  par_sel <- as.matrix(wcpo_pars$par_sel)
  double_normal_f <- data$sel_type_f == 2L
  par_sel[double_normal_f, 3:4] <- par_sel[double_normal_f, 3:4] - log(sd(data$len_mid))
  par_sel
}

make_parameters <- function(data, wcpo_pars, log_B0_start = log(cached_fit$report$B0),
                            par_sel = NULL, rdev_y = NULL) {
  if (is.null(par_sel)) par_sel <- make_par_sel_start(data, wcpo_pars)
  if (is.null(rdev_y)) rdev_y <- as.numeric(wcpo_pars$rdev_y)
  list(
    log_B0 = log_B0_start,
    log_h = as.numeric(wcpo_pars$log_h),
    log_sigma_r = as.numeric(wcpo_pars$log_sigma_r),
    log_cpue_q = rep(0, data$n_index),
    cpue_creep = rep(as.numeric(wcpo_pars$cpue_creep), data$n_index),
    log_cpue_tau = rep(log(cpue_tau), data$n_index),
    log_cpue_omega = rep(as.numeric(wcpo_pars$log_cpue_omega), data$n_index),
    log_lf_tau = as.numeric(log(rep(0.1, data$n_fishery))),
    log_wf_tau = rep(0, data$n_fishery),
    log_L1 = as.numeric(wcpo_pars$log_L1),
    log_L2 = as.numeric(wcpo_pars$log_L2),
    log_k = as.numeric(wcpo_pars$log_k),
    log_CV1 = as.numeric(wcpo_pars$log_CV1),
    log_CV2 = as.numeric(wcpo_pars$log_CV2),
    par_sel = par_sel,
    rdev_y = rdev_y
  )
}

make_current_map_sel <- function(parameters) {
  map_sel <- cached_fit$diagnostics$selectivity_raw$map_sel
  stopifnot("Cached selectivity map does not match par_sel dimensions" =
              identical(dim(map_sel), dim(parameters$par_sel)))
  map_sel
}

make_map <- function(data, parameters, estimate_selectivity = TRUE,
                     estimate_rdev = TRUE, fix_B0 = FALSE) {
  map_sel <- if (estimate_selectivity) {
    make_current_map_sel(parameters)
  } else {
    matrix(NA_integer_, nrow(parameters$par_sel), ncol(parameters$par_sel))
  }

  map <- list(
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
    par_sel = factor(map_sel),
    rdev_y = factor(if (estimate_rdev) seq_along(parameters$rdev_y) else rep(NA_integer_, length(parameters$rdev_y)))
  )
  if (fix_B0) map$log_B0 <- factor(NA)

  list(map = map, map_sel = map_sel)
}

make_bounds <- function(obj, data, map, map_sel) {
  lower <- rep(-Inf, length(obj$par))
  upper <- rep(Inf, length(obj$par))

  lower[names(obj$par) == "log_B0"] <- 12
  upper[names(obj$par) == "log_B0"] <- 22
  lower[names(obj$par) == "log_cpue_q"] <- log(0.1)
  upper[names(obj$par) == "log_cpue_q"] <- log(10)
  lower[names(obj$par) == "log_lf_tau"] <- -9
  upper[names(obj$par) == "log_lf_tau"] <- 9
  lower[names(obj$par) == "rdev_y"] <- -5
  upper[names(obj$par) == "rdev_y"] <- 5

  mu_len <- mean(data$len_mid)
  sd_len <- sd(data$len_mid)
  par_labels <- names(obj$par)
  par_sel_pos <- which(names(obj$par) == "par_sel")
  par_sel_level_pos <- setNames(par_sel_pos, levels(map$par_sel))

  for (f in seq_len(data$n_fishery)) {
    idx <- rep(NA_integer_, ncol(map_sel))
    for (p in seq_len(ncol(map_sel))) {
      level <- as.character(map_sel[f, p])
      if (!is.na(level)) {
        idx[p] <- par_sel_level_pos[[level]]
        par_labels[idx[p]] <- paste0("par_sel[", f, ",", p, "]")
      }
    }

    if (data$sel_type_f[f] == 1L) {
      if (!is.na(idx[1])) {
        lower[idx[1]] <- (5 - mu_len) / sd_len
        upper[idx[1]] <- (200 - mu_len) / sd_len
      }
      if (!is.na(idx[2])) {
        lower[idx[2]] <- log(4 / sd_len)
        upper[idx[2]] <- log(500 / sd_len)
      }
    } else {
      if (!is.na(idx[1])) {
        lower[idx[1]] <- (10.1 - mu_len) / sd_len
        upper[idx[1]] <- (200 - mu_len) / sd_len
      }
      if (!is.na(idx[2])) {
        lower[idx[2]] <- -7
        upper[idx[2]] <- 7
      }
      shift <- 2 * log(sd_len)
      if (!is.na(idx[3])) {
        lower[idx[3]] <- log(4) - shift
        upper[idx[3]] <- 8 - shift
      }
      if (!is.na(idx[4])) {
        lower[idx[4]] <- log(4) - shift
        upper[idx[4]] <- 8 - shift
      }
      if (!is.na(idx[5])) {
        lower[idx[5]] <- -9
        upper[idx[5]] <- 9
      }
      if (!is.na(idx[6])) {
        lower[idx[6]] <- -9
        upper[idx[6]] <- 9
      }
    }
  }

  list(lower = lower, upper = upper, labels = par_labels)
}

max_kkt_gradient <- function(par, grad, lower, upper, tol = 1e-7) {
  kkt <- abs(grad)
  at_lower <- is.finite(lower) & par <= lower + tol
  at_upper <- is.finite(upper) & par >= upper - tol
  kkt[at_lower] <- pmax(0, -grad[at_lower])
  kkt[at_upper] <- pmax(0, grad[at_upper])
  max(kkt)
}

prepare_data <- function(lf_switch = 1L, wf_switch = 1L,
                         lf_adjust_mult = 1, wf_adjust_mult = 1,
                         sex_ratio = NULL) {
  data <- wcpo_bet_data |>
    add_length_bins() |>
    set_quarter_cpue()

  lf_var_adjust <- rep(80 * lf_adjust_mult, data$n_fishery)
  lf_var_adjust[8] <- 5000 * lf_adjust_mult
  data <- prep_lf_data(
    data = data,
    lf_wide = lf_wide_global,
    lf_keep_fisheries = c(8, 9, 10, 11, 12, 13, 14),
    lf_var_adjust = lf_var_adjust,
    lf_switch = 1L
  )

  data$wt_bin_start <- 1
  data$wt_bin_width <- 1
  data$n_wt <- 200L
  wf_var_adjust <- rep(2000 * wf_adjust_mult, data$n_fishery)
  data <- prep_wf_data(
    data = data,
    wf_wide = wf_wide_global,
    wf_keep_fisheries = c(1, 2, 3, 4, 6, 7, 15),
    wf_switch = 1L,
    wf_var_adjust = wf_var_adjust
  )

  data$lf_switch <- as.integer(lf_switch)
  data$wf_switch <- as.integer(wf_switch)
  data$priors <- list()
  if (!is.null(sex_ratio)) data$sex_ratio <- rep(sex_ratio, data$n_age)

  data
}

baseline_start <- function(obj, remove_names = character()) {
  start <- cached_fit$opt$par
  if (length(remove_names)) {
    keep <- !names(start) %in% remove_names
    start <- start[keep]
  }
  if (length(start) == length(obj$par) && identical(names(start), names(obj$par))) {
    return(start)
  }
  obj$par
}

component_sum <- function(rep, name) {
  value <- rep[[name]]
  if (is.null(value)) return(NA_real_)
  sum(as.numeric(value))
}

summarise_fit <- function(label, scenario, obj, opt, bounds, data, rep, elapsed) {
  grad <- obj$gr(opt$par)
  q <- exp(opt$par[names(opt$par) == "log_cpue_q"])
  log_B0 <- log(as.numeric(rep$B0))
  lower_B0 <- bounds$lower[names(opt$par) == "log_B0"]
  upper_B0 <- bounds$upper[names(opt$par) == "log_B0"]
  if (!length(lower_B0)) lower_B0 <- NA_real_
  if (!length(upper_B0)) upper_B0 <- NA_real_

  tibble(
    label = label,
    scenario = scenario,
    convergence = opt$convergence,
    message = opt$message,
    objective = as.numeric(opt$objective),
    raw_max_gradient = max(abs(grad)),
    kkt_max_gradient = max_kkt_gradient(opt$par, grad, bounds$lower, bounds$upper),
    n_par = length(opt$par),
    elapsed_seconds = as.numeric(elapsed[["elapsed"]]),
    B0 = as.numeric(rep$B0),
    log_B0 = log_B0,
    log_B0_lower = lower_B0,
    log_B0_upper = upper_B0,
    final_spawning_biomass = tail(as.numeric(rep$spawning_biomass_y), 1),
    final_static_depletion = tail(as.numeric(rep$static_depletion_y), 1),
    final_dynamic_depletion = tail(as.numeric(rep$dynamic_depletion_y), 1),
    R0 = as.numeric(rep$R0),
    lp_prior = component_sum(rep, "lp_prior"),
    lp_penalty = component_sum(rep, "lp_penalty"),
    lp_rec = component_sum(rep, "lp_rec"),
    lp_cpue = component_sum(rep, "lp_cpue"),
    lp_lf = component_sum(rep, "lp_lf"),
    lp_wf = component_sum(rep, "lp_wf"),
    cpue_log_rmse = sqrt(mean((log(data$cpue_data$value) - log(as.numeric(rep$cpue_pred)))^2)),
    max_harvest_rate = max(as.numeric(rep$hrate_ysa), na.rm = TRUE),
    q_min = min(q),
    q_max = max(q),
    q_at_bound = any(abs(log(q) - log(0.1)) < 1e-5 | abs(log(q) - log(10)) < 1e-5),
    n_bounds = sum(
      (is.finite(bounds$lower) & abs(opt$par - bounds$lower) < 1e-5) |
        (is.finite(bounds$upper) & abs(opt$par - bounds$upper) < 1e-5)
    )
  )
}

run_fit <- function(label, scenario,
                    lf_switch = 1L, wf_switch = 1L,
                    lf_adjust_mult = 1, wf_adjust_mult = 1,
                    estimate_selectivity = TRUE,
                    fixed_selectivity = c("start", "retained"),
                    estimate_rdev = TRUE,
                    rdev_zero = FALSE,
                    fix_B0 = NULL,
                    sex_ratio = NULL,
                    max_restarts = 2L,
                    eval_max = 2000L,
                    iter_max = 2000L) {
  cat("\nSetting up ", label, ": ", scenario, "\n", sep = "")
  flush(stdout())
  fixed_selectivity <- match.arg(fixed_selectivity)
  data <- prepare_data(
    lf_switch = lf_switch,
    wf_switch = wf_switch,
    lf_adjust_mult = lf_adjust_mult,
    wf_adjust_mult = wf_adjust_mult,
    sex_ratio = sex_ratio
  )
  cat("  data prepared: lf_switch=", data$lf_switch,
      ", wf_switch=", data$wf_switch,
      ", n_lf=", length(data$lf_year),
      ", n_wf=", length(data$wf_year), "\n", sep = "")
  flush(stdout())

  par_sel <- if (fixed_selectivity == "retained") {
    cached_fit$diagnostics$selectivity_raw$par_est
  } else {
    NULL
  }
  rdev_y <- if (rdev_zero) rep(0, length(wcpo_bet_parameters$rdev_y)) else NULL
  log_B0_start <- if (is.null(fix_B0)) log(cached_fit$report$B0) else log(fix_B0)
  parameters <- make_parameters(
    data = data,
    wcpo_pars = wcpo_bet_parameters,
    log_B0_start = log_B0_start,
    par_sel = par_sel,
    rdev_y = rdev_y
  )
  map_info <- make_map(
    data = data,
    parameters = parameters,
    estimate_selectivity = estimate_selectivity,
    estimate_rdev = estimate_rdev,
    fix_B0 = !is.null(fix_B0)
  )
  cat("  building AD object...\n")
  flush(stdout())
  obj <- MakeADFun(func = cmb(opal_model, data), parameters = parameters, map = map_info$map)
  obj$env$tracemgc <- FALSE
  cat("  AD object built; setting bounds and start values...\n")
  flush(stdout())
  bounds <- make_bounds(obj, data, map_info$map, map_info$map_sel)

  remove_names <- character()
  if (!is.null(fix_B0)) remove_names <- c(remove_names, "log_B0")
  if (!estimate_selectivity) remove_names <- c(remove_names, "par_sel")
  if (!estimate_rdev) remove_names <- c(remove_names, "rdev_y")
  start <- baseline_start(obj, remove_names = remove_names)

  cat("Running ", label, ": start nll=", signif(obj$fn(start), 8),
      ", npar=", length(start), "\n", sep = "")
  flush(stdout())
  control <- list(eval.max = eval_max, iter.max = iter_max)
  timing <- system.time({
    opt <- list(par = start, convergence = NA_integer_, message = "not run", objective = obj$fn(start))
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
      kkt <- max_kkt_gradient(opt$par, grad, bounds$lower, bounds$upper)
      cat("  restart ", i, ": obj=", signif(opt$objective, 8),
          ", raw_grad=", signif(max(abs(grad)), 5),
          ", kkt_grad=", signif(kkt, 5),
          ", conv=", opt$convergence, "\n", sep = "")
      flush(stdout())
      if (kkt < 1e-3) break
    }
  })
  obj$env$last.par.best <- opt$par
  rep <- obj$report(opt$par)
  summary <- summarise_fit(label, scenario, obj, opt, bounds, data, rep, timing)

  list(
    label = label,
    scenario = scenario,
    summary = summary,
    opt = opt,
    report = rep,
    data = data,
    bounds = bounds,
    map_sel = map_info$map_sel
  )
}

xmodel_b0 <- xmodel_biomass |>
  mutate(B0 = ssb / depletion) |>
  group_by(model) |>
  summarise(B0 = median(B0, na.rm = TRUE), .groups = "drop")

ss3_fix_b0 <- xmodel_b0$B0[xmodel_b0$model == "02-fix-sel"]
mfcl_b0 <- xmodel_b0$B0[xmodel_b0$model == "MFCL-v11"]
ss3_base_b0 <- xmodel_b0$B0[xmodel_b0$model == "01-bet-base"]

cached_summary <- tibble(
  label = "cached_retained",
  scenario = "cached retained 15-fishery fit",
  convergence = cached_fit$opt$convergence,
  message = cached_fit$opt$message,
  objective = cached_fit$opt$objective,
  raw_max_gradient = NA_real_,
  kkt_max_gradient = cached_fit$kkt_grad %||% NA_real_,
  n_par = length(cached_fit$opt$par),
  elapsed_seconds = cached_fit$diagnostics$mle_runtime_seconds %||% NA_real_,
  B0 = as.numeric(cached_fit$report$B0),
  log_B0 = log(as.numeric(cached_fit$report$B0)),
  log_B0_lower = 12,
  log_B0_upper = 22,
  final_spawning_biomass = tail(as.numeric(cached_fit$report$spawning_biomass_y), 1),
  final_static_depletion = tail(as.numeric(cached_fit$report$static_depletion_y), 1),
  final_dynamic_depletion = tail(as.numeric(cached_fit$report$dynamic_depletion_y), 1),
  R0 = as.numeric(cached_fit$report$R0),
  lp_prior = component_sum(cached_fit$report, "lp_prior"),
  lp_penalty = component_sum(cached_fit$report, "lp_penalty"),
  lp_rec = component_sum(cached_fit$report, "lp_rec"),
  lp_cpue = component_sum(cached_fit$report, "lp_cpue"),
  lp_lf = component_sum(cached_fit$report, "lp_lf"),
  lp_wf = component_sum(cached_fit$report, "lp_wf"),
  cpue_log_rmse = sqrt(mean((log(cached_fit$data$cpue_data$value) - log(cached_fit$report$cpue_pred))^2)),
  max_harvest_rate = max(as.numeric(cached_fit$report$hrate_ysa), na.rm = TRUE),
  q_min = min(exp(cached_fit$opt$par[names(cached_fit$opt$par) == "log_cpue_q"])),
  q_max = max(exp(cached_fit$opt$par[names(cached_fit$opt$par) == "log_cpue_q"])),
  q_at_bound = FALSE,
  n_bounds = NA_integer_
)

scenarios <- list(
  list(
    label = "retained_refit",
    scenario = "current likelihood/map, start from cached optimum",
    estimate_selectivity = TRUE,
    fixed_selectivity = "start"
  ),
  list(
    label = "sel_fixed_retained",
    scenario = "current likelihood, selectivity fixed at retained estimates",
    estimate_selectivity = FALSE,
    fixed_selectivity = "retained"
  ),
  list(
    label = "sel_fixed_no_wf",
    scenario = "WF likelihood off, retained selectivity fixed",
    wf_switch = 0L,
    estimate_selectivity = FALSE,
    fixed_selectivity = "retained"
  ),
  list(
    label = "sel_fixed_no_lf",
    scenario = "LF likelihood off, retained selectivity fixed",
    lf_switch = 0L,
    estimate_selectivity = FALSE,
    fixed_selectivity = "retained"
  ),
  list(
    label = "sel_fixed_cpue_catch_only",
    scenario = "LF and WF likelihoods off, retained selectivity fixed",
    lf_switch = 0L,
    wf_switch = 0L,
    estimate_selectivity = FALSE,
    fixed_selectivity = "retained"
  ),
  list(
    label = "comp_downweight_10",
    scenario = "LF and WF effective sample sizes downweighted ten-fold",
    lf_adjust_mult = 10,
    wf_adjust_mult = 10,
    estimate_selectivity = TRUE,
    fixed_selectivity = "start"
  ),
  list(
    label = "no_rdev",
    scenario = "recruitment deviations fixed at zero",
    estimate_rdev = FALSE,
    rdev_zero = TRUE,
    estimate_selectivity = TRUE,
    fixed_selectivity = "start"
  ),
  list(
    label = "sex_ratio_half",
    scenario = "spawning-potential sex ratio set to 0.5",
    sex_ratio = 0.5,
    estimate_selectivity = TRUE,
    fixed_selectivity = "start"
  ),
  list(
    label = "profile_B0_SS3_02_fix_sel",
    scenario = "B0 fixed at SS3 02-fix-sel implied scale",
    fix_B0 = ss3_fix_b0,
    estimate_selectivity = TRUE,
    fixed_selectivity = "start"
  ),
  list(
    label = "profile_B0_MFCL",
    scenario = "B0 fixed at MFCL implied scale",
    fix_B0 = mfcl_b0,
    estimate_selectivity = TRUE,
    fixed_selectivity = "start"
  ),
  list(
    label = "profile_B0_SS3_base",
    scenario = "B0 fixed at SS3 base implied scale",
    fix_B0 = ss3_base_b0,
    estimate_selectivity = TRUE,
    fixed_selectivity = "start"
  )
)

custom_comp_mult <- Sys.getenv("OPAL_BET_COMP_DOWNWEIGHT_MULT", unset = "")
if (nzchar(custom_comp_mult)) {
  custom_comp_mult <- as.numeric(strsplit(custom_comp_mult, ",", fixed = TRUE)[[1]])
  custom_comp_mult <- custom_comp_mult[is.finite(custom_comp_mult) & custom_comp_mult > 0]
  scenarios <- c(
    scenarios,
    map(custom_comp_mult, function(mult) {
      label_mult <- gsub("\\.", "p", format(mult, trim = TRUE, scientific = FALSE))
      list(
        label = paste0("comp_downweight_", label_mult),
        scenario = paste0("LF and WF effective sample sizes downweighted ", mult, "-fold"),
        lf_adjust_mult = mult,
        wf_adjust_mult = mult,
        estimate_selectivity = TRUE,
        fixed_selectivity = "start"
      )
    })
  )
}

requested_labels <- Sys.getenv("OPAL_BET_SENS_LABELS", unset = "")
if (nzchar(requested_labels)) {
  keep_labels <- trimws(strsplit(requested_labels, ",", fixed = TRUE)[[1]])
  scenarios <- keep(map(scenarios, ~ if (.x$label %in% keep_labels) .x else NULL), Negate(is.null))
  if (!length(scenarios)) stop("No scenarios matched OPAL_BET_SENS_LABELS")
}

results <- vector("list", length(scenarios))
for (i in seq_along(scenarios)) {
  args <- scenarios[[i]]
  args$lf_switch <- args$lf_switch %||% 1L
  args$wf_switch <- args$wf_switch %||% 1L
  args$lf_adjust_mult <- args$lf_adjust_mult %||% 1
  args$wf_adjust_mult <- args$wf_adjust_mult %||% 1
  args$estimate_selectivity <- args$estimate_selectivity %||% TRUE
  args$fixed_selectivity <- args$fixed_selectivity %||% "start"
  args$estimate_rdev <- args$estimate_rdev %||% TRUE
  args$rdev_zero <- args$rdev_zero %||% FALSE
  args$fix_B0 <- args$fix_B0 %||% NULL
  args$sex_ratio <- args$sex_ratio %||% NULL
  results[[i]] <- do.call(run_fit, args)
}

summary_tbl <- bind_rows(cached_summary, map(results, "summary")) |>
  mutate(
    objective_delta_from_cached = objective - cached_fit$opt$objective,
    B0_ratio_to_cached = B0 / as.numeric(cached_fit$report$B0),
    B0_ratio_to_SS3_02_fix_sel = B0 / ss3_fix_b0,
    B0_ratio_to_MFCL = B0 / mfcl_b0
  )

out <- list(
  created = Sys.time(),
  opal_dir = opal_dir,
  xmodel_b0 = xmodel_b0,
  cached_summary = cached_summary,
  summary = summary_tbl,
  results = results
)

saveRDS(out, file.path(out_dir, "opal_bet_biomass_sensitivity.rds"))
write.csv(summary_tbl, file.path(out_dir, "opal_bet_biomass_sensitivity_summary.csv"), row.names = FALSE)

cat("\nSummary:\n")
print(summary_tbl |>
  select(label, objective, objective_delta_from_cached, B0, B0_ratio_to_cached,
         final_static_depletion, lp_penalty, lp_cpue, lp_lf, lp_wf,
         cpue_log_rmse, kkt_max_gradient, convergence))
