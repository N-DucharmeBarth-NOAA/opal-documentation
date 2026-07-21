#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(tibble)
})

msg <- function(...) {
  cat(format(Sys.time(), "%H:%M:%S"), ..., "\n")
  flush(stdout())
}

if (!requireNamespace("compResidual", quietly = TRUE)) {
  stop(
    "Package 'compResidual' is required to build the six-fishery residual cache. ",
    "Install it in a local maintenance environment; it is intentionally excluded from CI."
  )
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

file_provenance <- function(path) {
  info <- file.info(path)
  list(
    path = normalizePath(path, winslash = "/", mustWork = TRUE),
    size_bytes = unname(info$size),
    modified = as.character(info$mtime)
  )
}

raw_path <- Sys.getenv(
  "OPAL_BET_AGG6_DW0P1_RDS",
  "assets/cached/bet/agg6_dw0p1/agg6_result.rds"
)
summary_path <- Sys.getenv(
  "OPAL_BET_AGG6_SUMMARY_RDS",
  "assets/cached/bet/opal_bet_agg6_summary.rds"
)
out_path <- Sys.getenv(
  "OPAL_BET_AGG6_RESIDUALS_RDS",
  "assets/cached/bet/opal_bet_agg6_comp_residuals.rds"
)

for (path in c(raw_path, summary_path)) {
  if (!file.exists(path)) stop("Missing required source cache: ", path)
}

msg("Loading six-fishery raw result")
raw <- readRDS(raw_path)
if (is.null(raw$result$fit$data) || is.null(raw$result$report)) {
  stop("Raw cache does not contain result fit data and report: ", raw_path)
}
data <- raw$result$fit$data
report <- raw$result$report

msg("Loading six-fishery summary")
summary <- readRDS(summary_path)
agg6 <- summary$downweight0p1
if (is.null(agg6$grouping)) {
  stop("Summary cache does not contain downweight0p1 grouping: ", summary_path)
}
group_labels <- agg6$grouping %>%
  transmute(
    fishery_num = group,
    group_label = paste0("Group ", group, " (old F", old_fisheries, ")")
  )

ts_to_decimal_year <- function(ts) {
  data$first_yr + (as.integer(ts) - 1L) %/% data$n_season +
    ((as.integer(ts) - 1L) %% data$n_season) / data$n_season
}

composition_residual_rows <- function(type_label, obs_flat, pred_list,
                                      n_by_row, n_f, fishery_f,
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
      plot_year <- years[j]
      if (!is.finite(plot_year) || plot_year < 1000) {
        plot_year <- ts_to_decimal_year(years[j])
      }
      keep <- seq_along(residual)

      tibble(
        composition = type_label,
        type = paste(type_label, "composition"),
        Fishery = paste("Group", f),
        fishery_num = as.integer(f),
        row = as.integer(row_ids[j]),
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

msg("Computing length-composition residuals")
lf_residuals <- composition_residual_rows(
  type_label = "Length",
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
msg("Computing weight-composition residuals")
wf_residuals <- composition_residual_rows(
  type_label = "Weight",
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

rows <- bind_rows(lf_residuals, wf_residuals) %>%
  left_join(group_labels, by = "fishery_num")
required_columns <- c(
  "composition", "type", "Fishery", "fishery_num", "row", "year", "bin",
  "observed", "predicted", "n_eff", "residual", "sign", "group_label"
)
missing_columns <- setdiff(required_columns, names(rows))
if (length(missing_columns)) {
  stop("Residual rows are missing required columns: ", paste(missing_columns, collapse = ", "))
}
if (!nrow(rows) || anyNA(rows$fishery_num) || anyNA(rows$group_label)) {
  stop("Residual cache validation failed: empty rows or unmatched aggregation groups")
}

cache <- list(
  schema_version = 1L,
  residual_method = "compResidual::resMulti",
  created = Sys.time(),
  raw_cache = file_provenance(raw_path),
  summary_cache = file_provenance(summary_path),
  residual_rows = rows
)

dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
temp_path <- tempfile(pattern = "opal_bet_agg6_comp_residuals_", tmpdir = dirname(out_path), fileext = ".rds")
saveRDS(cache, temp_path)
if (!file.copy(temp_path, out_path, overwrite = TRUE)) {
  unlink(temp_path)
  stop("Failed to write residual cache: ", out_path)
}
unlink(temp_path)
msg("Wrote ", out_path, " (", nrow(rows), " residual cells)")
