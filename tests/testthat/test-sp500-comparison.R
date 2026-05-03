source_sp500_comparison <- function() {
  script <- system.file("examples", "sp500_compare_algorithms.R", package = "febama")
  if (!nzchar(script)) {
    script <- testthat::test_path(
      "..",
      "..",
      "inst",
      "examples",
      "sp500_compare_algorithms.R"
    )
  }
  env <- new.env(parent = globalenv())
  source(script, local = env)
  env
}

sp500_comparison_returns <- function(n = 42) {
  i <- seq_len(n)
  as.numeric(0.02 * sin(i / 3) + 0.01 * cos(i / 7) + 0.001 * i)
}

sp500_comparison_config <- function(env, alg_name = c("MAP", "sgld")) {
  alg_name <- match.arg(alg_name)
  env$sp500_paper_config(
    fore_model = c("naive_fore", "rw_drift_fore"),
    history_burn = 25,
    roll = NULL,
    feature_window = NULL,
    features_used = c("x_acf1", "entropy"),
    beta_shrinkage = 100,
    variable_selection = FALSE,
    init_optim = FALSE,
    alg_name = alg_name,
    n_iter = 2,
    max_batch_size = 4,
    n_epoch = 2,
    burnin_prop = 0,
    stepsize = 1e-4,
    lpd_features_parallel = FALSE,
    ncores = 1
  )
}

test_that("SP500 algorithm comparison returns MAP and SGLD rows", {
  sp500 <- source_sp500_comparison()

  returns <- sp500_comparison_returns()
  configs <- list(
    MAP = sp500_comparison_config(sp500, "MAP"),
    SGLD = sp500_comparison_config(sp500, "sgld")
  )

  result <- sp500$sp500_run_algorithm_comparison(
    returns = returns,
    n_origins = 1,
    configs = configs,
    seed = 10
  )
  summary <- sp500$sp500_algorithm_summary(result)
  delta <- sp500$sp500_algorithm_delta(summary)

  expect_equal(result$algorithm, c("MAP", "SGLD"))
  expect_named(summary, c("algorithm", "log_score", "mase", "smape", "elapsed_seconds"))
  expect_equal(delta$metric, c("log_score", "mase", "smape", "elapsed_seconds"))
  expect_true(all(is.finite(result$forecast)))
  expect_true(all(is.finite(result$log_score)))
  expect_true(all(is.finite(result$weight_naive_fore)))
  expect_true(all(is.finite(result$weight_rw_drift_fore)))
})

test_that("MAP outer updates keep all features when variable selection is disabled", {
  sp500 <- source_sp500_comparison()

  returns <- sp500_comparison_returns()
  config <- sp500_comparison_config(sp500, "MAP")
  series <- sp500$sp500_comparison_series(
    returns,
    origin = length(returns) - 1L,
    config = config
  )
  lpd_features <- clean_features(compute_lpd_features(series, config))
  fit <- fit_febama(lpd_features, config)

  expect_true(all(fit$betaIdx[[1]] == 1))
  expect_true(all(is.finite(fit$beta[[1]])))
})
