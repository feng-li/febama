# Compare MAP and SGLD inference on the S&P 500 example.
#
# Usage:
#   Rscript inst/examples/sp500_compare_algorithms.R data/sp500_daily_percent_log_returns.csv
#
# The script uses the same rolling-origin density/features for both algorithms
# at each origin, then reports one row per algorithm and origin. Set
# FEBAMA_SP500_FAST=1 to use the lightweight naive/random-walk forecasters.

sp500_comparison_script_path <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) > 0L) {
    return(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE))
  }

  for (frame in rev(sys.frames())) {
    if (!is.null(frame$ofile)) {
      return(normalizePath(frame$ofile, mustWork = FALSE))
    }
  }

  NA_character_
}

source_sp500_example <- function() {
  script <- system.file("examples", "sp500.R", package = "febama")
  if (!nzchar(script)) {
    this_script <- sp500_comparison_script_path()
    script <- if (is.na(this_script)) {
      file.path("inst", "examples", "sp500.R")
    } else {
      file.path(dirname(this_script), "sp500.R")
    }
  }

  source(script, local = parent.frame())
}

source_sp500_example()

numeric_env <- function(name, default, allow_na = FALSE) {
  value <- Sys.getenv(name, unset = as.character(default))
  if (!nzchar(value)) {
    return(default)
  }
  if (allow_na && tolower(value) == "na") {
    return(NA_real_)
  }

  out <- suppressWarnings(as.numeric(value))
  if (is.na(out)) {
    default
  } else {
    out
  }
}

sp500_validate_comparison_configs <- function(configs) {
  if (!is.list(configs) || length(configs) < 2L) {
    stop("`configs` must be a named list with at least two algorithm configs.",
         call. = FALSE)
  }
  if (is.null(names(configs)) || any(!nzchar(names(configs)))) {
    stop("`configs` must be named, for example list(MAP = ..., SGLD = ...).",
         call. = FALSE)
  }

  reference <- configs[[1L]]
  shared_fields <- c(
    "frequency",
    "ets_model",
    "ts_scale",
    "forecast_h",
    "train_h",
    "history_burn",
    "PI_level",
    "roll",
    "feature_window",
    "features_used",
    "fore_model"
  )

  for (i in seq_along(configs)) {
    config <- configs[[i]]
    if (config$algArgs$nIter < 2L) {
      stop("Algorithm comparison requires `n_iter >= 2` for every config.",
           call. = FALSE)
    }

    for (field in shared_fields) {
      if (!identical(config[[field]], reference[[field]])) {
        stop("Comparison configs must share `", field, "`.", call. = FALSE)
      }
    }
  }

  invisible(TRUE)
}

sp500_comparison_series <- function(returns, origin, config) {
  if (origin <= config$history_burn) {
    stop("`origin` must be larger than `config$history_burn`.", call. = FALSE)
  }
  if (origin + config$forecast_h > length(returns)) {
    stop("Not enough held-out returns for the requested forecast horizon.",
         call. = FALSE)
  }

  list(
    x = as.numeric(returns[seq_len(origin)]),
    xx = as.numeric(returns[seq.int(origin + 1L, origin + config$forecast_h)])
  )
}

sp500_missing_features <- function(lpd_features, config) {
  used <- unique(unlist(config$features_used, use.names = FALSE))
  used <- used[!is.na(used) & nzchar(used)]
  setdiff(used, colnames(lpd_features$feat))
}

sp500_compare_one_origin <- function(returns, origin, configs, seed = 1L) {
  sp500_validate_comparison_configs(configs)

  feature_config <- configs[[1L]]
  series <- sp500_comparison_series(returns, origin, feature_config)

  lpd_features <- compute_lpd_features(series, feature_config)
  lpd_features <- clean_features(lpd_features)
  missing_features <- sp500_missing_features(lpd_features, feature_config)
  if (length(missing_features) > 0L) {
    stop(
      "The configured S&P 500 feature set is missing from computed features: ",
      paste(missing_features, collapse = ", "),
      call. = FALSE
    )
  }

  rows <- lapply(seq_along(configs), function(i) {
    algorithm <- names(configs)[[i]]
    config <- configs[[i]]

    if (!is.null(seed)) {
      set.seed(seed + origin)
    }
    started <- proc.time()[["elapsed"]]
    fit <- fit_febama(lpd_features, config)
    elapsed <- proc.time()[["elapsed"]] - started

    forecast <- forecast_febama(
      data = series,
      config = config,
      lpd_features = lpd_features,
      fit = fit
    )

    weights <- as.numeric(forecast$w_time_varying[, 1L])
    names(weights) <- paste0("weight_", rownames(forecast$w_time_varying))

    cbind(
      data.frame(
        algorithm = algorithm,
        origin = origin,
        forecast = as.numeric(forecast$ff_feature),
        actual = as.numeric(forecast$xx),
        log_score = as.numeric(forecast$err_feature[, "lpds"]),
        mase = as.numeric(forecast$err_feature[, "mase_err_h"]),
        smape = as.numeric(forecast$err_feature[, "smape_err_h"]),
        elapsed_seconds = as.numeric(elapsed),
        check.names = FALSE
      ),
      as.data.frame(as.list(weights), check.names = FALSE)
    )
  })

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

sp500_run_algorithm_comparison <- function(returns,
                                           n_origins = 1L,
                                           final_origin = length(returns) - 1L,
                                           configs,
                                           seed = 1L) {
  sp500_validate_comparison_configs(configs)
  if (configs[[1L]]$forecast_h != 1L) {
    stop("The S&P 500 paper experiment uses one-step forecasts; set forecast_h = 1.",
         call. = FALSE)
  }

  origins <- seq.int(final_origin - n_origins + 1L, final_origin)
  rows <- lapply(
    origins,
    sp500_compare_one_origin,
    returns = returns,
    configs = configs,
    seed = seed
  )

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

sp500_algorithm_summary <- function(results) {
  stats::aggregate(
    cbind(log_score, mase, smape, elapsed_seconds) ~ algorithm,
    data = results,
    FUN = mean
  )
}

sp500_algorithm_delta <- function(summary,
                                  candidate = "SGLD",
                                  reference = "MAP") {
  if (!all(c(candidate, reference) %in% summary$algorithm)) {
    return(data.frame())
  }

  metrics <- setdiff(names(summary), "algorithm")
  candidate_row <- summary[summary$algorithm == candidate, metrics, drop = FALSE]
  reference_row <- summary[summary$algorithm == reference, metrics, drop = FALSE]
  data.frame(
    metric = metrics,
    candidate = candidate,
    reference = reference,
    difference = as.numeric(candidate_row[1L, ] - reference_row[1L, ]),
    check.names = FALSE
  )
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  path <- if (length(args) > 0L) args[[1L]] else Sys.getenv("FEBAMA_SP500_CSV")
  if (!nzchar(path)) {
    stop(
      "Provide a CSV path: Rscript inst/examples/sp500_compare_algorithms.R path/to/sp500.csv",
      call. = FALSE
    )
  }

  fast_mode <- truthy_env("FEBAMA_SP500_FAST")
  fore_model <- if (fast_mode) {
    c("naive_fore", "rw_drift_fore")
  } else {
    c("garch_fore", "rgarch_fore", "sv_fore")
  }
  if (!fast_mode) {
    sp500_check_optional_packages(fore_model)
  }

  forecast_h <- integer_env("FEBAMA_SP500_FORECAST_H", 1L)
  n_origins <- integer_env("FEBAMA_SP500_ORIGINS", 1L)
  ncores <- integer_env("FEBAMA_SP500_NCORES", max(1L, parallel::detectCores() - 1L))
  variable_selection <- truthy_env("FEBAMA_SP500_VARIABLE_SELECTION")

  default_n_iter <- integer_env("FEBAMA_SP500_N_ITER", 2L)
  map_n_iter <- integer_env("FEBAMA_SP500_MAP_N_ITER", default_n_iter)
  sgld_n_iter <- integer_env("FEBAMA_SP500_SGLD_N_ITER", default_n_iter)
  sgld_max_batch_size <- integer_env("FEBAMA_SP500_SGLD_MAX_BATCH_SIZE", 108L)
  sgld_n_epoch <- integer_env("FEBAMA_SP500_SGLD_N_EPOCH", 10L)
  sgld_burnin_prop <- numeric_env("FEBAMA_SP500_SGLD_BURNIN_PROP", 0.4)
  sgld_stepsize <- numeric_env("FEBAMA_SP500_SGLD_STEPSIZE", 0.1, allow_na = TRUE)
  sgld_gamma <- numeric_env("FEBAMA_SP500_SGLD_GAMMA", 0.55)
  sgld_a <- numeric_env("FEBAMA_SP500_SGLD_A", 0.4)
  sgld_b <- numeric_env("FEBAMA_SP500_SGLD_B", 10)
  seed <- integer_env("FEBAMA_SP500_SEED", 1L)

  common_args <- list(
    fore_model = fore_model,
    forecast_h = forecast_h,
    variable_selection = variable_selection,
    lpd_features_parallel = ncores > 1L,
    ncores = ncores
  )

  configs <- list(
    MAP = do.call(
      sp500_paper_config,
      c(common_args, list(alg_name = "MAP", n_iter = map_n_iter))
    ),
    SGLD = do.call(
      sp500_paper_config,
      c(
        common_args,
        list(
          alg_name = "sgld",
          n_iter = sgld_n_iter,
          max_batch_size = sgld_max_batch_size,
          n_epoch = sgld_n_epoch,
          burnin_prop = sgld_burnin_prop,
          stepsize = sgld_stepsize,
          sgld_gamma = sgld_gamma,
          sgld_a = sgld_a,
          sgld_b = sgld_b
        )
      )
    )
  )

  returns <- sp500_read_returns(path)
  results <- sp500_run_algorithm_comparison(
    returns = returns,
    n_origins = n_origins,
    configs = configs,
    seed = seed
  )
  summary <- sp500_algorithm_summary(results)
  delta <- sp500_algorithm_delta(summary)

  cat("\nPer-origin MAP vs SGLD results:\n")
  print(results)
  cat("\nAverage by algorithm:\n")
  print(summary)
  if (nrow(delta) > 0L) {
    cat("\nSGLD minus MAP:\n")
    print(delta)
  }

  output <- if (length(args) > 1L) {
    args[[2L]]
  } else {
    Sys.getenv("FEBAMA_SP500_COMPARISON_CSV")
  }
  if (nzchar(output)) {
    dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(results, output, row.names = FALSE)
    message("Wrote per-origin comparison to ", output)
  }

  invisible(list(results = results, summary = summary, delta = delta))
}

if (sys.nframe() == 0L) {
  main()
}
