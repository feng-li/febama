# Paper-style S&P 500 daily return example for FEBAMA.
#
# Usage:
#   Rscript inst/examples/sp500.R path/to/sp500.csv
#
# The CSV must contain either:
#   - a daily percent log-return column named return, returns, log_return,
#     log_returns, ret, or sp500_return; or
#   - a price column named adjusted, adj_close, adj.close, close, price, or value.
#
# The default settings mirror Section 4 of docs/paper.pdf:
#   - one-step forecasts;
#   - 1,250 trading-day rolling windows for volatility models;
#   - 100-day sliding windows for feature calculation;
#   - GARCH, realized GARCH, and stochastic-volatility base models;
#   - the 15 stock-market features listed in Table 3.
#
# Full out-of-sample replication is computationally expensive because each
# origin refits the historical density/feature training set. Set
# FEBAMA_SP500_ORIGINS to run more than the default final-origin example.

example_script_path <- function() {
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

load_febama_for_example <- function() {
  if (requireNamespace("febama", quietly = TRUE)) {
    library(febama)
    return(invisible(TRUE))
  }

  script <- example_script_path()
  repo_root <- if (is.na(script)) {
    normalizePath(".", mustWork = FALSE)
  } else {
    normalizePath(file.path(dirname(script), "..", ".."), mustWork = FALSE)
  }

  if (file.exists(file.path(repo_root, "DESCRIPTION")) &&
      requireNamespace("pkgload", quietly = TRUE)) {
    getExportedValue("pkgload", "load_all")(repo_root, quiet = TRUE)
    return(invisible(TRUE))
  }

  stop(
    "Install febama first, or install pkgload to run this example from a source checkout.",
    call. = FALSE
  )
}

load_febama_for_example()

sp500_feature_set <- function() {
  c(
    "alpha",
    "arch_acf",
    "arch_r2",
    "beta",
    "crossing_points",
    "diff1x_pacf5",
    "diff2_acf1",
    "diff2_acf10",
    "entropy",
    "garch_acf",
    "garch_r2",
    "nonlinearity",
    "trend",
    "unitroot_kpss",
    "x_acf1"
  )
}

sp500_paper_config <- function(fore_model = c("garch_fore", "rgarch_fore", "sv_fore"),
                               forecast_h = 1,
                               history_burn = 1250,
                               roll = 1250,
                               feature_window = 100,
                               features_used = sp500_feature_set(),
                               beta_shrinkage = 1000,
                               variable_selection = FALSE,
                               varsel_init = "all-in",
                               init_optim = TRUE,
                               alg_name = c("MAP", "sgld"),
                               n_iter = 1,
                               max_batch_size = 108,
                               n_epoch = 10,
                               burnin_prop = 0.4,
                               stepsize = 0.1,
                               sgld_gamma = 0.55,
                               sgld_a = 0.4,
                               sgld_b = 10,
                               lpd_features_parallel = TRUE,
                               ncores = max(1L, parallel::detectCores() - 1L)) {
  alg_name <- match.arg(alg_name)

  febama_config(
    frequency = 1,
    forecast_h = forecast_h,
    train_h = 1,
    history_burn = history_burn,
    roll = roll,
    feature_window = feature_window,
    PI_level = 90,
    fore_model = fore_model,
    features_used = features_used,
    lpd_features_parallel = lpd_features_parallel,
    ncores = ncores,
    variable_selection = variable_selection,
    varsel_init = varsel_init,
    beta_shrinkage = beta_shrinkage,
    init_optim = init_optim,
    alg_name = alg_name,
    n_iter = n_iter,
    max_batch_size = max_batch_size,
    n_epoch = n_epoch,
    burnin_prop = burnin_prop,
    stepsize = stepsize,
    sgld_gamma = sgld_gamma,
    sgld_a = sgld_a,
    sgld_b = sgld_b
  )
}

sp500_read_returns <- function(path) {
  if (!file.exists(path)) {
    stop("Cannot find S&P 500 input file: ", path, call. = FALSE)
  }

  data <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  original_names <- names(data)
  normalized_names <- normalize_column_names(original_names)

  date_col <- first_matching_column(normalized_names, c("date", "time", "timestamp"))
  if (!is.na(date_col)) {
    data <- data[order(as.Date(data[[date_col]])), , drop = FALSE]
  }

  return_col <- first_matching_column(
    normalized_names,
    c("return", "returns", "log_return", "log_returns", "ret", "sp500_return")
  )
  if (!is.na(return_col)) {
    returns <- as.numeric(data[[return_col]])
    return(stats::na.omit(returns))
  }

  price_col <- first_matching_column(
    normalized_names,
    c("adjusted", "adjusted_close", "adj_close", "adjclose", "adj.close",
      "close", "price", "value")
  )
  if (is.na(price_col)) {
    stop(
      "Input CSV must contain a return column or a close/adjusted price column.",
      call. = FALSE
    )
  }

  prices <- as.numeric(data[[price_col]])
  returns <- 100 * diff(log(prices))
  stats::na.omit(returns)
}

sp500_run_one_step <- function(returns,
                               origin = length(returns) - 1L,
                               config = sp500_paper_config()) {
  if (origin <= config$history_burn) {
    stop("`origin` must be larger than `config$history_burn`.", call. = FALSE)
  }
  if (origin + config$forecast_h > length(returns)) {
    stop("Not enough held-out returns for the requested forecast horizon.", call. = FALSE)
  }

  series <- list(
    x = as.numeric(returns[seq_len(origin)]),
    xx = as.numeric(returns[seq.int(origin + 1L, origin + config$forecast_h)])
  )

  lpd_features <- compute_lpd_features(series, config)
  lpd_features <- clean_features(lpd_features)
  missing_features <- setdiff(config$features_used[[1]], colnames(lpd_features$feat))
  if (length(missing_features) > 0L) {
    stop(
      "The configured S&P 500 feature set is missing from computed features: ",
      paste(missing_features, collapse = ", "),
      call. = FALSE
    )
  }

  fit <- fit_febama(lpd_features, config)
  forecast <- forecast_febama(
    data = series,
    config = config,
    lpd_features = lpd_features,
    fit = fit
  )

  list(
    origin = origin,
    config = config,
    lpd_features = lpd_features,
    fit = fit,
    weights = prepare_febama_weights(fit, config, burnin = 0),
    forecast = forecast,
    performance = summarize_performance(forecast)
  )
}

sp500_run_backtest <- function(returns,
                               n_origins = 1L,
                               final_origin = length(returns) - 1L,
                               config = sp500_paper_config()) {
  if (config$forecast_h != 1L) {
    stop("The S&P 500 paper experiment uses one-step forecasts; set forecast_h = 1.",
         call. = FALSE)
  }

  origins <- seq.int(final_origin - n_origins + 1L, final_origin)
  rows <- lapply(origins, function(origin) {
    result <- sp500_run_one_step(returns, origin = origin, config = config)
    data.frame(
      origin = origin,
      forecast = as.numeric(result$forecast$ff_feature),
      actual = as.numeric(result$forecast$xx),
      log_score = as.numeric(result$forecast$err_feature[, "lpds"]),
      mase = as.numeric(result$forecast$err_feature[, "mase_err_h"]),
      smape = as.numeric(result$forecast$err_feature[, "smape_err_h"]),
      t(result$forecast$w_time_varying[, 1, drop = FALSE]),
      check.names = FALSE
    )
  })

  do.call(rbind, rows)
}

sp500_check_optional_packages <- function(fore_model) {
  package_map <- list(
    garch_fore = "rugarch",
    egarch_fore = "rugarch",
    rgarch_fore = c("rugarch", "highfrequency", "xts"),
    sv_fore = "stochvol"
  )
  required <- unique(unlist(package_map[intersect(names(package_map), fore_model)]))
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0L) {
    stop(
      "The paper-style S&P 500 models require optional packages: ",
      paste(missing, collapse = ", "),
      ". Install them or set FEBAMA_SP500_FAST=1 for a plumbing-only example.",
      call. = FALSE
    )
  }
}

normalize_column_names <- function(x) {
  tolower(gsub("[^[:alnum:]]+", "_", x))
}

first_matching_column <- function(normalized_names, candidates) {
  idx <- match(candidates, normalized_names, nomatch = NA_integer_)
  idx <- idx[!is.na(idx)]
  if (length(idx) == 0L) {
    return(NA_integer_)
  }
  idx[[1]]
}

truthy_env <- function(name, default = FALSE) {
  value <- Sys.getenv(name, unset = if (default) "1" else "0")
  tolower(value) %in% c("1", "true", "yes", "y")
}

integer_env <- function(name, default) {
  value <- Sys.getenv(name, unset = as.character(default))
  out <- suppressWarnings(as.integer(value))
  if (is.na(out)) {
    default
  } else {
    out
  }
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  path <- if (length(args) > 0L) args[[1L]] else Sys.getenv("FEBAMA_SP500_CSV")
  if (!nzchar(path)) {
    stop(
      "Provide a CSV path: Rscript inst/examples/sp500.R path/to/sp500.csv",
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
  n_iter <- integer_env("FEBAMA_SP500_N_ITER", 1L)
  variable_selection <- truthy_env("FEBAMA_SP500_VARIABLE_SELECTION")

  config <- sp500_paper_config(
    fore_model = fore_model,
    forecast_h = forecast_h,
    variable_selection = variable_selection,
    n_iter = n_iter,
    lpd_features_parallel = ncores > 1L,
    ncores = ncores
  )
  returns <- sp500_read_returns(path)
  result <- sp500_run_backtest(
    returns = returns,
    n_origins = n_origins,
    config = config
  )

  print(result)
  invisible(result)
}

if (sys.nframe() == 0L) {
  main()
}
