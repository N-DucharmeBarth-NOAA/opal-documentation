#!/usr/bin/env Rscript

cat(format(Sys.time(), "%H:%M:%S"), "Starting BET cache promotion\n")
flush(stdout())

extra_r_libs <- Sys.getenv("OPAL_EXTRA_R_LIBS", unset = "")
extra_r_lib_paths <- if (nzchar(extra_r_libs)) {
  strsplit(extra_r_libs, .Platform$path.sep, fixed = TRUE)[[1]]
} else {
  character()
}
r_minor <- strsplit(R.version$minor, ".", fixed = TRUE)[[1]][1]
default_user_lib <- file.path(
  path.expand("~"),
  "R",
  paste0(R.version$platform, "-library"),
  paste(R.version$major, r_minor, sep = ".")
)
extra_r_lib_paths <- unique(c(extra_r_lib_paths, default_user_lib))
extra_r_lib_paths <- extra_r_lib_paths[dir.exists(extra_r_lib_paths)]
if (length(extra_r_lib_paths)) .libPaths(c(extra_r_lib_paths, .libPaths()))

suppressPackageStartupMessages({
  library(RTMB)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
})

if (!requireNamespace("compResidual", quietly = TRUE)) {
  stop("Package 'compResidual' is required to rebuild BET composition residuals")
}

msg <- function(...) {
  cat(format(Sys.time(), "%H:%M:%S"), ..., "\n")
  flush(stdout())
}

normalise_composition <- function(x) {
  total <- sum(x, na.rm = TRUE)
  if (!is.finite(total) || total <= 0) return(rep(NA_real_, length(x)))
  x / total
}

quiet_comp_residual_warnings <- function(expr) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      if (grepl("NA/NaN function evaluation", conditionMessage(w), fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

sdnr_summary <- function(x, conf = 0.95) {
  x <- x[is.finite(x)]
  df <- length(x) - 1L
  est <- stats::sd(x)
  if (df < 1L) {
    return(tibble(SDNR = est, LCI = NA_real_, HCI = NA_real_))
  }
  alpha <- (1 - conf) / 2
  tibble(
    SDNR = est,
    LCI = sqrt(df * est^2 / stats::qchisq(1 - alpha, df)),
    HCI = sqrt(df * est^2 / stats::qchisq(alpha, df))
  )
}

repo_dir <- Sys.getenv("OPAL_DOCUMENTATION_DIR", "/home/darcy/Projects/opal-documentation")
opal_dir <- Sys.getenv("OPAL_SOURCE_DIR", "/home/darcy/Projects/opal")
scenario <- Sys.getenv("OPAL_BET_PROMOTE_SCENARIO", "directional_moderate")
cache_dir <- file.path(repo_dir, "assets", "cached", "bet")
candidate_path <- file.path(cache_dir, paste0("bet_reweight_", scenario), "wf_index_results.rds")
out_path <- file.path(cache_dir, "opal_bet_fit.rds")

if (!file.exists(candidate_path)) stop("Missing candidate: ", candidate_path)

for (file in list.files(file.path(opal_dir, "R"), pattern = "\\.R$", full.names = TRUE)) {
  sys.source(file, envir = .GlobalEnv)
}
load(file.path(opal_dir, "data", "wcpo_bet_lf.rda"))
load(file.path(opal_dir, "data", "wcpo_bet_wf.rda"))

candidate <- readRDS(candidate_path)[[1]]
obj <- candidate$fit$obj
data <- candidate$fit$data
data$wcpo_bet_lf <- wcpo_bet_lf
data$wcpo_bet_wf <- wcpo_bet_wf
obj$env$last.par.best <- candidate$opt$par
obj$env$last.par <- candidate$opt$par
report <- obj$report(candidate$opt$par)
par_list <- obj$env$parList(candidate$opt$par)

msg("Loaded candidate ", scenario, "; objective=", signif(candidate$opt$objective, 8))

make_component_obj <- function(cpue_switch, lf_switch, wf_switch) {
  data_i <- data
  data_i$cpue_switch <- as.integer(cpue_switch)
  data_i$lf_switch <- as.integer(lf_switch)
  data_i$wf_switch <- as.integer(wf_switch)
  obj_i <- MakeADFun(
    func = cmb(opal_model, data_i),
    parameters = candidate$fit$parameters,
    map = candidate$fit$map
  )
  obj_i$env$tracemgc <- FALSE
  obj_i
}

msg("Computing objective decomposition")
obj_none <- make_component_obj(0L, 0L, 0L)
obj_cpue <- make_component_obj(1L, 0L, 0L)
obj_lf <- make_component_obj(0L, 1L, 0L)
obj_wf <- make_component_obj(0L, 0L, 1L)
none_nll <- obj_none$fn(candidate$opt$par)
none_report <- obj_none$report(candidate$opt$par)
objective_components <- c(
  "Prior contribution" = sum(as.numeric(none_report$lp_prior)),
  "Recruitment-deviation contribution" =
    sum(as.numeric(none_report$lp_rec)) + sum(as.numeric(none_report$lp_init_rec)),
  "CPUE contribution" = obj_cpue$fn(candidate$opt$par) - none_nll,
  "Length-composition contribution" = obj_lf$fn(candidate$opt$par) - none_nll,
  "Weight-composition contribution" = obj_wf$fn(candidate$opt$par) - none_nll,
  "Harvest penalty" = sum(as.numeric(none_report$lp_penalty))
)
component_error <- sum(objective_components) - candidate$opt$objective
if (!is.finite(component_error) || abs(component_error) > 1e-5) {
  stop("Objective components do not sum to optimizer objective; difference = ", component_error)
}

msg("Simulating CPUE")
set.seed(9010)
n_cpue_sim <- 20L
cpue_sim <- map_dfr(seq_len(n_cpue_sim), function(s) {
  sim <- obj$simulate(candidate$opt$par)
  tibble(
    ts = data$cpue_data$ts,
    year = data$cpue_data$year,
    month = data$cpue_data$month,
    index = data$cpue_data$index,
    sim_id = s,
    simulated_log = as.numeric(sim$cpue_log_obs)
  )
})
report$cpue_sim <- cpue_sim

msg("Computing CPUE OSA residuals")
osa_cpue <- oneStepPredict(
  obj = obj,
  observation.name = "cpue_log_obs",
  method = "oneStepGeneric",
  trace = FALSE
)
report$osa_cpue_residual <- as.numeric(osa_cpue$res)
cpue_sdnr <- bind_rows(
  tibble(index = "All indices", residual = report$osa_cpue_residual),
  tibble(index = month.name[data$cpue_data$month], residual = report$osa_cpue_residual)
) |>
  group_by(.data$index) |>
  group_modify(~ bind_cols(
    tibble(
      n = sum(is.finite(.x$residual)),
      MAR = median(abs(.x$residual), na.rm = TRUE)
    ),
    sdnr_summary(.x$residual)
  )) |>
  ungroup() |>
  mutate(index = factor(
    .data$index,
    levels = c("All indices", month.name[sort(unique(data$cpue_data$month))])
  )) |>
  arrange(.data$index) |>
  mutate(index = as.character(.data$index))

flat_to_count_list <- function(v, n_f, fishery_f, minbin, maxbin) {
  out <- vector("list", length(n_f))
  offset <- 0L
  for (k in seq_along(n_f)) {
    f <- fishery_f[k]
    nbins <- maxbin[f] - minbin[f] + 1L
    n_values <- n_f[k] * nbins
    out[[k]] <- matrix(
      v[seq.int(offset + 1L, length.out = n_values)],
      nrow = n_f[k],
      ncol = nbins,
      byrow = TRUE
    )
    offset <- offset + n_values
  }
  out
}

ts_to_decimal_year <- function(ts) {
  data$first_yr + (as.integer(ts) - 1L) %/% data$n_season +
    ((as.integer(ts) - 1L) %% data$n_season) / data$n_season
}

composition_residual_rows <- function(type_label, y_label, obs_flat,
                                      pred_list, n_by_row, n_f, fishery_f,
                                      year_fi, minbin, maxbin, bins) {
  obs_list <- flat_to_count_list(obs_flat, n_f, fishery_f, minbin, maxbin)
  row_offset <- 0L

  map_dfr(seq_along(pred_list), function(i) {
    f <- fishery_f[i]
    n_rows <- n_f[i]
    row_ids <- seq.int(row_offset + 1L, length.out = n_rows)
    row_offset <<- row_offset + n_rows
    active_bins <- bins[minbin[f]:maxbin[f]]
    pred <- pred_list[[i]]
    obs <- obs_list[[i]]
    years <- year_fi[[as.character(f)]]

    msg("  ", type_label, " fishery ", f, " residuals")
    map_dfr(seq_len(n_rows), function(j) {
      obs_prop <- normalise_composition(obs[j, ])
      pred_prop <- normalise_composition(pred[j, ])
      n_eff <- n_by_row[row_ids[j]]

      if (anyNA(obs_prop) || anyNA(pred_prop) ||
          !is.finite(n_eff) || n_eff <= 0) {
        return(tibble())
      }

      pred_prop <- pred_prop + 1e-6
      pred_prop <- pred_prop / sum(pred_prop)
      residual <- as.numeric(quiet_comp_residual_warnings(
        compResidual::resMulti(
          obs = matrix(obs_prop * n_eff, ncol = 1L),
          pred = matrix(pred_prop, ncol = 1L)
        )
      ))
      keep <- seq_along(residual)
      plot_year <- years[j]
      if (!is.finite(plot_year) || plot_year < 1000) {
        plot_year <- ts_to_decimal_year(years[j])
      }

      tibble(
        composition = type_label,
        type = paste(type_label, "composition"),
        Fishery = paste("Fishery", f),
        fishery_num = f,
        row = row_ids[j],
        year = plot_year,
        bin = active_bins[keep],
        observed = obs_prop[keep],
        predicted = pred_prop[keep],
        n_eff = n_eff,
        residual = residual,
        sign = if_else(residual >= 0, "Positive", "Negative")
      )
    })
  })
}

msg("Computing composition OSA residuals")
lf_residuals <- composition_residual_rows(
  type_label = "Length",
  y_label = "Length (cm)",
  obs_flat = data$lf_obs_flat,
  pred_list = report$lf_pred,
  n_by_row = data$lf_n,
  n_f = data$lf_n_f,
  fishery_f = data$lf_fishery_f,
  year_fi = data$lf_year_fi,
  minbin = data$lf_minbin,
  maxbin = data$lf_maxbin,
  bins = data$len_mid
)
wf_residuals <- composition_residual_rows(
  type_label = "Weight",
  y_label = "Weight (kg)",
  obs_flat = data$wf_obs_flat,
  pred_list = report$wf_pred,
  n_by_row = data$wf_n,
  n_f = data$wf_n_f,
  fishery_f = data$wf_fishery_f,
  year_fi = data$wf_year_fi,
  minbin = data$wf_minbin,
  maxbin = data$wf_maxbin,
  bins = data$wt_mid
)

composition_sdnr <- bind_rows(lf_residuals, wf_residuals) |>
  filter(is.finite(.data$residual)) |>
  group_by(.data$composition, .data$Fishery, .data$fishery_num) |>
  group_modify(~ bind_cols(
    tibble(
      n = nrow(.x),
      MAR = median(abs(.x$residual), na.rm = TRUE)
    ),
    sdnr_summary(.x$residual)
  )) |>
  ungroup()

simulate_row_mean <- function(pred_list, fishery_f, row_fishery, n_eff, nsim = 20L) {
  map(seq_len(nsim), function(s) {
    map(seq_along(pred_list), function(i) {
      f <- fishery_f[i]
      pred <- as.matrix(pred_list[[i]])
      n_vec <- as.numeric(n_eff[row_fishery == f])
      sim_rows <- t(vapply(seq_len(nrow(pred)), function(r) {
        p <- as.numeric(pred[r, ])
        p[!is.finite(p) | p < 0] <- 0
        p <- if (sum(p) > 0) p / sum(p) else rep(1 / length(p), length(p))
        size <- max(1L, as.integer(round(n_vec[r])))
        as.numeric(rmultinom(1L, size = size, prob = p)) / size
      }, numeric(ncol(pred))))
      colMeans(sim_rows, na.rm = TRUE)
    })
  })
}
set.seed(9020)
report$lf_sim_mean <- simulate_row_mean(report$lf_pred, data$lf_fishery_f, data$lf_fishery, data$lf_n)
report$wf_sim_mean <- simulate_row_mean(report$wf_pred, data$wf_fishery_f, data$wf_fishery, data$wf_n)

msg("Computing biomass uncertainty")
hess_ndeps <- 1e-4
hess <- stats::optimHess(
  par = candidate$opt$par,
  fn = obj$fn,
  gr = obj$gr,
  control = list(ndeps = rep(hess_ndeps, length(candidate$opt$par)))
)
hess_eigen <- eigen((hess + t(hess)) / 2, symmetric = TRUE, only.values = TRUE)$values
sdr <- tryCatch(
  RTMB::sdreport(obj, hessian.fixed = hess),
  error = function(e) {
    msg("sdreport rejected supplied Hessian: ", conditionMessage(e))
    RTMB::sdreport(obj)
  }
)
sdr_report <- summary(sdr, "report")
get_report_se <- function(name, n) {
  rows <- which(rownames(sdr_report) == name)
  if (!length(rows)) {
    return(tibble(estimate = rep(NA_real_, n), se = rep(NA_real_, n)))
  }
  rows <- rows[seq_len(min(length(rows), n))]
  tibble(
    estimate = sdr_report[rows, "Estimate"],
    se = sdr_report[rows, "Std. Error"]
  )
}
n_bio <- length(report$spawning_biomass_y)
sb <- get_report_se("spawning_biomass_y", n_bio)
sb0 <- get_report_se("spawning_biomass0_y", n_bio)
static_dep <- get_report_se("static_depletion_y", n_bio)
dyn_dep <- get_report_se("dynamic_depletion_y", n_bio)
bio_year <- data$cpue_data$year[1] + (seq_len(n_bio) - 1) / 4
biomass_uncertainty <- tibble(
  year = bio_year,
  spawning_biomass = as.numeric(report$spawning_biomass_y),
  spawning_biomass_se = sb$se,
  sb_sb0 = as.numeric(report$static_depletion_y),
  sb_sb0_se = static_dep$se,
  dynamic_depletion = as.numeric(report$dynamic_depletion_y),
  dynamic_depletion_se = dyn_dep$se,
  spawning_biomass0 = as.numeric(report$spawning_biomass0_y),
  spawning_biomass0_se = sb0$se
)

msg("Building selectivity and catch diagnostics")
catch_table <- as.data.frame.table(data$catch_obs_ysf, responseName = "catch")
names(catch_table)[seq_len(3)] <- c("ts", "season", "fishery")
catch_by_fishery <- catch_table |>
  transmute(
    ts = as.integer(.data$ts),
    season = as.integer(.data$season),
    fishery = as.integer(.data$fishery),
    catch = as.numeric(.data$catch),
    predicted_catch = as.numeric(report$catch_pred_ysf[cbind(
      as.integer(.data$ts),
      as.integer(.data$season),
      as.integer(.data$fishery)
    )])
  ) |>
  left_join(
    data$cpue_data |> distinct(.data$ts, .data$year, .data$month),
    by = "ts"
  ) |>
  mutate(
    decimal_year = .data$year + (.data$month - 1) / 12,
    units_code = data$catch_units_f[.data$fishery],
    units = if_else(.data$units_code == 1, "weight", "numbers")
  )

diagnostics <- list(
  cpue_sdnr = cpue_sdnr,
  composition_sdnr = composition_sdnr,
  lf_residuals = lf_residuals,
  wf_residuals = wf_residuals,
  residual_source = "`compResidual::resMulti`",
  objective_components = objective_components,
  objective_component_total = sum(objective_components),
  objective_component_full_objective = candidate$opt$objective,
  objective_component_method = "AD objective switch increments at the fitted parameter vector",
  composition_simulation_method = paste(
    "Fitted-model row-level multinomial simulations;",
    "sub-unit effective sample sizes rounded up for display before row means are averaged"
  ),
  mle_runtime_seconds = unname(candidate$runtime[["elapsed"]]),
  mle_runtime_platform = paste(R.version$platform, R.version$version.string),
  mle_runtime_control = list(
    log_B0_lower = 12,
    hessian_ndeps = hess_ndeps,
    cpue_tau = if (!is.null(candidate$cpue_tau)) {
      candidate$cpue_tau
    } else {
      unique(round(sqrt(pmax(
        as.numeric(report$cpue_sigma)^2 - data$cpue_data$se^2, 0
      )), 10))
    },
    sel_specs = candidate$sel_specs,
    lf_var_adjust_scalar = candidate$lf_var_adjust_scalar,
    lf_var_adjust_override = candidate$lf_var_adjust_override,
    wf_var_adjust = candidate$wf_var_adjust,
    wf_var_adjust_override = candidate$wf_var_adjust_override
  ),
  mle_runtime_iterations = candidate$opt$iterations,
  mle_runtime_evaluations = candidate$opt$evaluations,
  mle_runtime_objective = candidate$opt$objective,
  mle_runtime_convergence = candidate$opt$convergence,
  biomass_uncertainty = biomass_uncertainty,
  biomass_uncertainty_method =
    "RTMB::sdreport with supplied optimHess fixed-effect Hessian, ndeps = 1e-4",
  biomass_uncertainty_pdHess = isTRUE(sdr$pdHess),
  biomass_uncertainty_min_hessian_eigenvalue = min(hess_eigen, na.rm = TRUE),
  recruitment_raw = list(
    rdev_y = as.numeric(par_list$rdev_y),
    map_rdev = candidate$fit$map$rdev_y,
    cpue_lookup = data$cpue_data |> select(.data$ts, .data$year, .data$month)
  ),
  selectivity_raw = list(
    par_est = par_list$par_sel,
    par_ref = candidate$fit$parameters$par_sel,
    sel_type_f = data$sel_type_f,
    len_mid = data$len_mid,
    sel_specs = candidate$sel_specs,
    map_sel = candidate$fit$map_sel,
    opt_par = candidate$opt$par,
    bounds = candidate$fit$bounds[c("lower", "upper")]
  ),
  catch_by_fishery = catch_by_fishery
)

promoted <- list(
  source = candidate_path,
  label = candidate$label,
  opt = candidate$opt,
  raw_grad = candidate$raw_grad,
  kkt_grad = candidate$kkt_grad,
  data = data,
  report = report,
  diagnostics = diagnostics
)

saveRDS(promoted, out_path)
msg("Wrote ", out_path)
msg("B0=", round(report$B0), ", final SB/SB0=", signif(tail(report$spawning_biomass_y, 1) / report$B0, 4),
    ", pdHess=", isTRUE(sdr$pdHess), ", min Hessian eigen=", signif(min(hess_eigen), 4))
