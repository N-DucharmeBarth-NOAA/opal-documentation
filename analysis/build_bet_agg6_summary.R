#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
})

msg <- function(...) {
  cat(format(Sys.time(), "%H:%M:%S"), ..., "\n")
  flush(stdout())
}

sel_logistic_plain <- function(len, par) {
  mu <- mean(len)
  sd_len <- sd(len)
  a50 <- mu + par[1] * sd_len
  w95 <- exp(par[2]) * sd_len
  1 / (1 + exp(-log(19) * (len - a50) / w95))
}

sel_double_normal_plain <- function(len, par) {
  mu <- mean(len)
  sd_len <- sd(len)
  peak <- mu + par[1] * sd_len
  upselex <- exp(par[3]) * sd_len^2
  downselex <- exp(par[4]) * sd_len^2
  point1 <- 1 / (1 + exp(-par[5]))
  point2 <- 1 / (1 + exp(-par[6]))
  j2 <- length(len)
  bin_width <- len[2] - len[1]
  peak2 <- peak + bin_width + (0.99 * len[j2] - peak - bin_width) /
    (1 + exp(-par[2]))
  t1min <- exp(-(len[1] - peak)^2 / upselex)
  t2min <- exp(-(len[j2] - peak2)^2 / downselex)
  t1 <- len - peak
  t2 <- len - peak2
  join1 <- 1 / (1 + exp(-(20 / (1 + abs(t1))) * t1))
  join2 <- 1 / (1 + exp(-(20 / (1 + abs(t2))) * t2))
  asc <- point1 + (1 - point1) * (exp(-t1^2 / upselex) - t1min) /
    (1 - t1min)
  dsc <- 1 + (point2 - 1) * (exp(-t2^2 / downselex) - 1) /
    (t2min - 1)
  asc * (1 - join1) + join1 * (1 - join2 + dsc * join2)
}

composition_mean <- function(data, report, type = c("lf", "wf")) {
  type <- match.arg(type)
  if (type == "lf") {
    fishery <- data$lf_fishery
    fishery_f <- data$lf_fishery_f
    bins <- data$len_mid
    obs <- data$lf_obs_in
    pred <- report$lf_pred
    bin_name <- "length"
    x_label <- "Length (cm)"
  } else {
    fishery <- data$wf_fishery
    fishery_f <- data$wf_fishery_f
    bins <- data$wt_mid
    obs <- data$wf_obs_in
    pred <- report$wf_pred
    bin_name <- "weight"
    x_label <- "Weight (kg)"
  }

  obs_df <- as_tibble(obs, .name_repair = ~ as.character(bins)) %>%
    mutate(row = row_number(), fishery = fishery) %>%
    pivot_longer(
      cols = -c(row, fishery),
      names_to = bin_name,
      values_to = "proportion"
    ) %>%
    mutate(
      "{bin_name}" := as.numeric(.data[[bin_name]]),
      series = "Observed"
    )

  pred_df <- map_dfr(seq_along(pred), function(i) {
    as_tibble(pred[[i]], .name_repair = ~ as.character(bins)) %>%
      mutate(row = row_number(), fishery = fishery_f[i]) %>%
      pivot_longer(
        cols = -c(row, fishery),
        names_to = bin_name,
        values_to = "proportion"
      ) %>%
      mutate(
        "{bin_name}" := as.numeric(.data[[bin_name]]),
        series = "Predicted"
      )
  })

  bind_rows(obs_df, pred_df) %>%
    group_by(fishery, .data[[bin_name]], series) %>%
    summarise(proportion = mean(proportion, na.rm = TRUE), .groups = "drop") %>%
    mutate(
      composition = if_else(type == "lf", "Length", "Weight"),
      x_label = x_label
    )
}

selectivity_curves <- function(data, par_sel) {
  map_dfr(seq_len(data$n_fishery), function(f) {
    sel <- if (data$sel_type_f[f] == 1L) {
      sel_logistic_plain(data$len_mid, par_sel[f, ])
    } else {
      sel_double_normal_plain(data$len_mid, par_sel[f, ])
    }
    tibble(
      fishery = f,
      length = data$len_mid,
      selectivity = as.numeric(sel),
      form = if_else(data$sel_type_f[f] == 1L, "logistic", "double-normal")
    )
  })
}

catch_by_group <- function(data) {
  ts_lookup <- data$cpue_data %>%
    distinct(ts, year, month) %>%
    mutate(decimal_year = year + (month - 1) / 12)
  as.data.frame.table(data$catch_obs_ysf, responseName = "catch") %>%
    transmute(
      ts = as.integer(Var1),
      season = as.integer(Var2),
      fishery = as.integer(Var3),
      catch = as.numeric(catch)
    ) %>%
    left_join(ts_lookup, by = "ts") %>%
    mutate(
      units_code = data$catch_units_f[fishery],
      units = if_else(units_code == 1, "weight", "numbers")
    )
}

load_run <- function(path) {
  msg("loading ", path)
  x <- readRDS(path)
  result <- x$result
  data <- result$fit$data
  report <- result$report
  msg("building CPUE and biomass summaries")
  cpue <- data$cpue_data %>%
    mutate(
      index_month = factor(
        month,
        levels = sort(unique(month)),
        labels = month.name[sort(unique(month))]
      ),
      decimal_year = year + (month - 1) / 12,
      pred = as.numeric(report$cpue_pred),
      sigma = as.numeric(report$cpue_sigma),
      lower = exp(log(value) - sigma),
      upper = exp(log(value) + sigma)
    )
  biomass <- tibble(
    ts = seq_along(report$spawning_biomass_y),
    year = data$cpue_data$year[1] + (ts - 1) / 4,
    spawning_biomass = as.numeric(report$spawning_biomass_y),
    spawning_biomass0 = as.numeric(report$spawning_biomass0_y),
    sb_sb0 = spawning_biomass / as.numeric(report$B0),
    dynamic_depletion = as.numeric(report$dynamic_depletion_y)
  )
  lf_adjust_override <- result$lf_var_adjust_override
  msg("building composition, catch, and selectivity summaries")
  list(
    grouping = x$grouping,
    old_to_new = x$old_to_new,
    metrics = x$metrics,
    lf_metrics = x$lf_metrics,
    wf_metrics = x$wf_metrics,
    cpue = cpue,
    biomass = biomass,
    catch = catch_by_group(data),
    lf_mean = composition_mean(data, report, "lf"),
    wf_mean = composition_mean(data, report, "wf"),
    selectivity = selectivity_curves(data, result$par_list$par_sel),
    sel_status = x$metrics$sel_status,
    weighting = tibble(
      cpue_tau = if (!is.null(result$cpue_tau)) {
        result$cpue_tau
      } else {
        unique(round(sqrt(pmax(
          as.numeric(report$cpue_sigma)^2 - data$cpue_data$se^2, 0
        )), 10))[1]
      },
      wf_var_adjust = result$wf_var_adjust,
      lf_group1_var_adjust = if (length(lf_adjust_override)) {
        unname(lf_adjust_override[[1]])
      } else {
        NA_real_
      }
    )
  )
}

paths <- c(
  full_g5 = Sys.getenv("OPAL_BET_AGG6_BASE_RDS", "/tmp/opal_bet_agg6/agg6_result.rds"),
  no_desc = Sys.getenv("OPAL_BET_AGG6_NO_DESC_RDS", "/tmp/opal_bet_agg6_g5_no_desc/agg6_result.rds"),
  downweight0p1 = Sys.getenv(
    "OPAL_BET_AGG6_DW0P1_RDS",
    "assets/cached/bet/agg6_dw0p1/agg6_result.rds"
  ),
  downweight0p25 = Sys.getenv(
    "OPAL_BET_AGG6_DW0P25_RDS",
    "assets/cached/bet/agg6_dw0p25/agg6_result.rds"
  ),
  downweight0p5 = Sys.getenv(
    "OPAL_BET_AGG6_DW0P5_RDS",
    "assets/cached/bet/agg6_dw0p5/agg6_result.rds"
  ),
  upweight2 = Sys.getenv(
    "OPAL_BET_AGG6_UPW2_RDS",
    "assets/cached/bet/agg6_upw2/agg6_result.rds"
  ),
  upweight4 = Sys.getenv(
    "OPAL_BET_AGG6_UPW4_RDS",
    "assets/cached/bet/agg6_upw4/agg6_result.rds"
  )
)

if (!file.exists(paths[["full_g5"]])) stop("Missing ", paths[["full_g5"]])

msg("building aggregated six-fishery summary")
summary <- list(
  full_g5 = load_run(paths[["full_g5"]]),
  no_desc = if (file.exists(paths[["no_desc"]])) load_run(paths[["no_desc"]]) else NULL,
  downweight0p1 = if (file.exists(paths[["downweight0p1"]])) {
    load_run(paths[["downweight0p1"]])
  } else {
    NULL
  },
  downweight0p25 = if (file.exists(paths[["downweight0p25"]])) {
    load_run(paths[["downweight0p25"]])
  } else {
    NULL
  },
  downweight0p5 = if (file.exists(paths[["downweight0p5"]])) {
    load_run(paths[["downweight0p5"]])
  } else {
    NULL
  },
  upweight2 = if (file.exists(paths[["upweight2"]])) load_run(paths[["upweight2"]]) else NULL,
  upweight4 = if (file.exists(paths[["upweight4"]])) load_run(paths[["upweight4"]]) else NULL,
  created = Sys.time()
)

out_path <- Sys.getenv(
  "OPAL_BET_AGG6_SUMMARY_RDS",
  "assets/cached/bet/opal_bet_agg6_summary.rds"
)
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
saveRDS(summary, out_path)
cat("Wrote", out_path, "\n")
