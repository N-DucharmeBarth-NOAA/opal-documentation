#!/usr/bin/env Rscript

cat(format(Sys.time(), "%H:%M:%S"), "Starting R package load\n")
flush(stdout())

suppressPackageStartupMessages({
  library(RTMB)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(ggplot2)
})

msg <- function(...) {
  cat(format(Sys.time(), "%H:%M:%S"), ..., "\n")
  flush(stdout())
}

parse_scenario_names <- function(value, default) {
  if (!nzchar(value)) return(default)
  strsplit(value, ",", fixed = TRUE)[[1]]
}

fmt_named <- function(x) {
  if (!length(x)) return("")
  paste(paste(names(x), x, sep = ":"), collapse = ";")
}

repo_dir <- Sys.getenv("OPAL_DOCUMENTATION_DIR", "/home/darcy/Projects/opal-documentation")
opal_dir <- Sys.getenv("OPAL_SOURCE_DIR", "/home/darcy/Projects/opal")
cache_dir <- file.path(repo_dir, "assets", "cached", "bet")
cpue_tau <- as.numeric(Sys.getenv("OPAL_BET_CPUE_TAU", "1e-8"))
if (!is.finite(cpue_tau) || cpue_tau <= 0) {
  stop("OPAL_BET_CPUE_TAU must be a positive finite value")
}
msg("Starting BET composition reweight candidate setup")
base_fit_path <- file.path(cache_dir, "opal_bet_fit.rds")
base_fit <- readRDS(base_fit_path)

if (!nzchar(Sys.getenv("OPAL_BET_LF_LOG_B0_LOWER", unset = ""))) {
  Sys.setenv(OPAL_BET_LF_LOG_B0_LOWER = "12")
}

scenarios <- list(
  directional_moderate = list(
    lf_var_adjust_scalar = 400,
    lf_var_adjust_override = c(`8` = 25000, `10` = 1600),
    wf_var_adjust_scalar = 10000,
    wf_var_adjust_override = c(`3` = 40000, `6` = 40000, `7` = 40000)
  ),
  directional_mid = list(
    lf_var_adjust_scalar = 500,
    lf_var_adjust_override = c(`8` = 30000, `10` = 2000),
    wf_var_adjust_scalar = 8000,
    wf_var_adjust_override = c(`3` = 60000, `6` = 60000, `7` = 60000)
  ),
  directional_strong = list(
    lf_var_adjust_scalar = 300,
    lf_var_adjust_override = c(`8` = 20000, `10` = 2400),
    wf_var_adjust_scalar = 5000,
    wf_var_adjust_override = c(`3` = 80000, `6` = 80000, `7` = 80000)
  )
)

scenario_names <- parse_scenario_names(
  Sys.getenv("OPAL_BET_REWEIGHT_SCENARIOS", unset = ""),
  names(scenarios)
)
unknown <- setdiff(scenario_names, names(scenarios))
if (length(unknown)) {
  stop("Unknown scenario(s): ", paste(unknown, collapse = ", "))
}

setwd(opal_dir)
msg("Sourcing opal BET diagnostic helpers")
sys.source("inst/scripts/bet_wf_index_diagnostics.R", globalenv())
combo$helper$source_opal_r()
inputs <- combo$helper$load_bet_inputs()
load("data/wcpo_bet_wf.rda")
msg("Loaded opal BET inputs")

wf_wide <- make_wf_wide(wcpo_bet_wf)
lf_wide <- combo$helper$make_lf_wide(inputs$lf)
sel_specs <- base_fit$diagnostics$selectivity_raw$sel_specs
lf_fisheries <- c(8L, 9L, 10L, 11L, 12L, 13L, 14L)
wf_fisheries <- c(1L, 2L, 3L, 4L, 6L, 7L, 15L)

prepare_targeted_lf_wf_data <- function(data, lf_wide, wf_wide,
                                        lf_var_adjust_scalar,
                                        lf_var_adjust_override,
                                        wf_var_adjust_scalar,
                                        wf_var_adjust_override) {
  data <- combo$helper$add_length_bins(data)
  data <- combo$helper$set_quarter_cpue(data)

  lf_var_adjust <- rep(lf_var_adjust_scalar, data$n_fishery)
  if (length(lf_var_adjust_override)) {
    lf_var_adjust[as.integer(names(lf_var_adjust_override))] <-
      as.numeric(lf_var_adjust_override)
  }
  data <- prep_lf_data(
    data = data,
    lf_wide = lf_wide,
    lf_keep_fisheries = lf_fisheries,
    lf_var_adjust = lf_var_adjust,
    lf_switch = 1L
  )

  data$wt_bin_start <- 1
  data$wt_bin_width <- 1
  data$n_wt <- 200L
  wf_var_adjust <- rep(wf_var_adjust_scalar, data$n_fishery)
  if (length(wf_var_adjust_override)) {
    wf_var_adjust[as.integer(names(wf_var_adjust_override))] <-
      as.numeric(wf_var_adjust_override)
  }
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
  data
}

build_targeted_fit <- function(scenario) {
  data <- prepare_targeted_lf_wf_data(
    data = inputs$data,
    lf_wide = lf_wide,
    wf_wide = wf_wide,
    lf_var_adjust_scalar = scenario$lf_var_adjust_scalar,
    lf_var_adjust_override = scenario$lf_var_adjust_override,
    wf_var_adjust_scalar = scenario$wf_var_adjust_scalar,
    wf_var_adjust_override = scenario$wf_var_adjust_override
  )
  parameters <- combo$helper$make_parameters(data, inputs$parameters)
  parameters$log_cpue_tau <- rep(log(cpue_tau), data$n_index)
  map_info <- make_index_map(
    data = data,
    parameters = parameters,
    estimate_index_sel = TRUE,
    index_sel_cols = c(1L),
    sel_specs = sel_specs
  )
  obj <- MakeADFun(
    func = cmb(opal_model, data),
    parameters = parameters,
    map = map_info$map
  )
  obj$env$tracemgc <- FALSE
  bounds <- combo$helper$make_bounds(obj, data, map_info$map, map_info$map_sel)

  list(
    data = data,
    parameters = parameters,
    map = map_info$map,
    map_sel = map_info$map_sel,
    obj = obj,
    bounds = bounds
  )
}

simulate_row_mean_composition <- function(pred_list, fishery_f, row_fishery,
                                          bins, n_eff, nsim = 20L) {
  map_dfr(seq_along(pred_list), function(i) {
    fishery <- fishery_f[i]
    pred <- as.matrix(pred_list[[i]])
    n_vec <- as.numeric(n_eff[row_fishery == fishery])
    if (length(n_vec) != nrow(pred)) {
      stop("Effective-N rows do not match predicted rows for fishery ", fishery)
    }

    map_dfr(seq_len(nsim), function(s) {
      sim_rows <- t(vapply(seq_len(nrow(pred)), function(r) {
        p <- as.numeric(pred[r, ])
        p[!is.finite(p) | p < 0] <- 0
        p <- if (sum(p) > 0) p / sum(p) else rep(1 / length(p), length(p))
        size <- max(1L, as.integer(round(n_vec[r])))
        as.numeric(rmultinom(1L, size = size, prob = p)) / size
      }, numeric(ncol(pred))))

      tibble(
        fishery = paste("Fishery", fishery),
        bin = bins,
        proportion = colMeans(sim_rows, na.rm = TRUE),
        sim_id = s
      )
    })
  })
}

make_mean_composition <- function(data, report, kind = c("lf", "wf")) {
  kind <- match.arg(kind)
  if (kind == "lf") {
    pred <- report$lf_pred
    fishery_f <- data$lf_fishery_f
    row_fishery <- data$lf_fishery
    bins <- data$len_mid
    obs_matrix <- data$lf_obs_in
  } else {
    pred <- report$wf_pred
    fishery_f <- data$wf_fishery_f
    row_fishery <- data$wf_fishery
    bins <- data$wt_mid
    obs_matrix <- data$wf_obs_in
  }

  pred_df <- map_dfr(seq_along(pred), function(i) {
    as_tibble(pred[[i]], .name_repair = ~ as.character(bins)) |>
      mutate(row = row_number(), fishery = paste("Fishery", fishery_f[i])) |>
      pivot_longer(
        cols = -c(row, fishery),
        names_to = "bin",
        values_to = "proportion"
      ) |>
      mutate(bin = as.numeric(.data$bin), series = "Predicted")
  })

  obs_df <- as_tibble(obs_matrix, .name_repair = ~ as.character(bins)) |>
    mutate(row = row_number(), fishery = paste("Fishery", row_fishery)) |>
    pivot_longer(
      cols = -c(row, fishery),
      names_to = "bin",
      values_to = "proportion"
    ) |>
    mutate(bin = as.numeric(.data$bin), series = "Observed")

  bind_rows(obs_df, pred_df) |>
    group_by(.data$fishery, .data$bin, .data$series) |>
    summarise(proportion = mean(.data$proportion, na.rm = TRUE), .groups = "drop")
}

plot_composition_preview <- function(result, out_path) {
  data <- result$fit$data
  report <- result$report
  role <- c(Observed = "#b22222", Predicted = "#d7301f", Simulated = "#1f78b4")

  one_plot <- function(kind) {
    if (kind == "lf") {
      pred <- report$lf_pred
      fishery_f <- data$lf_fishery_f
      row_fishery <- data$lf_fishery
      bins <- data$len_mid
      n_eff <- data$lf_n
      xlab <- "Length (cm)"
      title <- "Length composition"
    } else {
      pred <- report$wf_pred
      fishery_f <- data$wf_fishery_f
      row_fishery <- data$wf_fishery
      bins <- data$wt_mid
      n_eff <- data$wf_n
      xlab <- "Weight (kg)"
      title <- "Weight composition"
    }

    mean_df <- make_mean_composition(data, report, kind)
    levels <- paste("Fishery", sort(unique(fishery_f)))
    mean_df <- mean_df |>
      mutate(fishery = factor(.data$fishery, levels = levels))
    sim_df <- simulate_row_mean_composition(
      pred_list = pred,
      fishery_f = fishery_f,
      row_fishery = row_fishery,
      bins = bins,
      n_eff = n_eff,
      nsim = 20L
    ) |>
      mutate(fishery = factor(.data$fishery, levels = levels))
    total_neff <- tibble(
      fishery = paste("Fishery", row_fishery),
      effective_n = as.numeric(n_eff)
    ) |>
      group_by(.data$fishery) |>
      summarise(effective_n = sum(.data$effective_n, na.rm = TRUE), .groups = "drop") |>
      mutate(
        fishery = factor(.data$fishery, levels = levels),
        label = sprintf("Total Neff = %.1f", .data$effective_n)
      )

    ggplot() +
      geom_col(
        data = filter(mean_df, .data$series == "Observed"),
        aes(x = .data$bin, y = .data$proportion),
        fill = role[["Observed"]],
        width = 1.8
      ) +
      geom_line(
        data = sim_df,
        aes(x = .data$bin, y = .data$proportion, group = .data$sim_id),
        colour = role[["Simulated"]],
        alpha = 0.35,
        linewidth = 0.35
      ) +
      geom_line(
        data = filter(mean_df, .data$series == "Predicted"),
        aes(x = .data$bin, y = .data$proportion),
        colour = role[["Predicted"]],
        linewidth = 0.8
      ) +
      geom_text(
        data = total_neff,
        aes(x = Inf, y = Inf, label = .data$label),
        hjust = 1.05,
        vjust = 1.25,
        size = 2.5,
        inherit.aes = FALSE
      ) +
      facet_wrap(~ fishery, ncol = 2, scales = "free_y") +
      scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
      labs(title = title, x = xlab, y = "Mean proportion") +
      theme_bw() +
      theme(panel.grid.minor = element_blank())
  }

  pdf(out_path, width = 8, height = 7)
  print(one_plot("lf"))
  print(one_plot("wf"))
  dev.off()
}

fit_scenario <- function(name, scenario) {
  out_dir <- file.path(cache_dir, paste0("bet_reweight_", name))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  raw_path <- file.path(out_dir, "wf_index_results.rds")
  msg("Running ", name)
  msg(
    "LF scalar=", scenario$lf_var_adjust_scalar,
    ", LF overrides=", fmt_named(scenario$lf_var_adjust_override),
    ", WF scalar=", scenario$wf_var_adjust_scalar,
    ", WF overrides=", fmt_named(scenario$wf_var_adjust_override)
  )

  if (identical(Sys.getenv("OPAL_BET_REWEIGHT_USE_CACHE", unset = ""), "true") &&
      file.exists(raw_path)) {
    msg("Using cached candidate result from ", raw_path)
    result <- readRDS(raw_path)[[1]]
    plot_composition_preview(result, file.path(out_dir, "composition_preview.pdf"))
    return(result)
  }

  built <- build_targeted_fit(scenario)
  obj <- built$obj
  start <- obj$par
  if (length(base_fit$opt$par) == length(start) &&
      identical(names(base_fit$opt$par), names(start))) {
    start <- base_fit$opt$par
    msg("Warm-starting from retained BET cache")
  }

  control <- list(
    eval.max = as.integer(Sys.getenv("OPAL_BET_REWEIGHT_EVAL_MAX", "2000")),
    iter.max = as.integer(Sys.getenv("OPAL_BET_REWEIGHT_ITER_MAX", "2000"))
  )
  max_restarts <- as.integer(Sys.getenv("OPAL_BET_REWEIGHT_MAX_RESTARTS", "4"))
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
  par_list <- obj$env$parList(opt$par)
  grad <- obj$gr(opt$par)
  kkt_grad <- combo$max_kkt_gradient(
    opt$par, grad, built$bounds$lower, built$bounds$upper
  )

  result <- list(
    label = name,
    lf_var_adjust_scalar = scenario$lf_var_adjust_scalar,
    wf_var_adjust = scenario$wf_var_adjust_scalar,
    lf_var_adjust_override = scenario$lf_var_adjust_override,
    wf_var_adjust_override = scenario$wf_var_adjust_override,
    estimate_index_sel = TRUE,
    index_sel_cols = c(1L),
    lf_fisheries = lf_fisheries,
    wf_fisheries = wf_fisheries,
    cpue_tau = cpue_tau,
    sel_specs = sel_specs,
    fit = built,
    opt = opt,
    report = report,
    par_list = par_list,
    raw_grad = max(abs(grad)),
    kkt_grad = kkt_grad,
    runtime = runtime
  )

  saveRDS(list(result), raw_path)
  metrics <- make_metrics(result) |>
    mutate(
      lf_var_adjust_scalar = scenario$lf_var_adjust_scalar,
      runtime_seconds = unname(runtime[["elapsed"]])
    )
  write.csv(metrics, file.path(out_dir, "wf_index_metrics.csv"), row.names = FALSE)
  write.csv(make_lf_metrics(result), file.path(out_dir, "wf_index_lf_metrics.csv"), row.names = FALSE)
  write.csv(make_wf_metrics(result), file.path(out_dir, "wf_index_wf_metrics.csv"), row.names = FALSE)
  plot_composition_preview(result, file.path(out_dir, "composition_preview.pdf"))

  msg("Wrote ", out_dir)
  print(metrics |>
    select(
      label, n_par, n_lf_obs, n_wf_obs, convergence, nll,
      raw_max_gradient, kkt_max_gradient, B0, lp_cpue, lp_lf, lp_wf,
      log_cpue_rmse, final_spawning_biomass_ratio, runtime_seconds
    ))
  invisible(result)
}

results <- imap(scenarios[scenario_names], function(scenario, name) {
  fit_scenario(name, scenario)
})
saveRDS(results, file.path(cache_dir, "bet_reweight_candidates.rds"))
msg("Wrote ", file.path(cache_dir, "bet_reweight_candidates.rds"))
