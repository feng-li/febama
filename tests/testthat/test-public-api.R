tiny_series <- function(n = 36) {
  i <- seq_len(n)
  as.numeric(0.05 * i + sin(i / 3) + cos(i / 7))
}

tiny_config <- function(lpd_features_parallel = FALSE, ncores = 1) {
  febama_config(
    frequency = 1,
    forecast_h = 2,
    train_h = 1,
    history_burn = 25,
    fore_model = c("naive_fore", "rw_drift_fore"),
    features_used = c("x_acf1", "entropy"),
    lpd_features_parallel = lpd_features_parallel,
    ncores = ncores,
    variable_selection = FALSE,
    init_optim = FALSE,
    n_iter = 1
  )
}

test_that("febama_config validates and normalizes public options", {
  config <- febama_config(
    frequency = 1,
    forecast_h = 2,
    train_h = 1,
    history_burn = 25,
    fore_model = c("naive_fore", "rw_drift_fore", "auto.arima_fore"),
    features_used = list("x_acf1", c("entropy", "unitroot_kpss")),
    lpd_features_parallel = TRUE,
    ncores = 2,
    variable_selection = FALSE,
    init_optim = FALSE
  )

  expect_equal(config$frequency, 1)
  expect_equal(config$forecast_h, 2)
  expect_equal(config$lpd_features_parl, list(par = TRUE, ncores = 2))
  expect_length(config$features_used, 2)
  expect_equal(config$features_used[[1]], "x_acf1")
  expect_equal(config$features_used[[2]], c("entropy", "unitroot_kpss"))
  expect_null(config$varSelArgs[[1]]$cand)

  expect_error(
    febama_config(fore_model = "naive_fore"),
    "at least two forecasting functions"
  )
  expect_error(
    febama_config(
      fore_model = c("naive_fore", "rw_drift_fore", "auto.arima_fore"),
      features_used = list("x_acf1")
    ),
    "one entry per non-baseline model"
  )
})

test_that("clean_features removes NA feature columns and keeps scaling metadata", {
  raw_features <- cbind(
    keep_a = c(1, 2, 3),
    drop_b = c(1, NA, 3),
    keep_c = c(3, 4, 5)
  )
  scaled_features <- scale(raw_features)
  lpd_features <- list(
    lpd = matrix(log(0.5), nrow = 3, ncol = 2),
    feat = scaled_features
  )

  cleaned <- clean_features(lpd_features)

  expect_equal(colnames(cleaned$feat), c("keep_a", "keep_c"))
  expect_named(cleaned$feat_mean, c("keep_a", "keep_c"))
  expect_named(cleaned$feat_sd, c("keep_a", "keep_c"))
  expect_equal(unname(cleaned$feat_mean), c(2, 4))
  expect_equal(unname(cleaned$feat_sd), c(1, 1))
  expect_false(anyNA(cleaned$feat))
})

test_that("public fitting and forecasting wrappers return expected shapes", {
  set.seed(1)
  values <- tiny_series(38)
  series <- list(x = values[1:36], xx = values[37:38])
  config <- tiny_config()

  lpd_features <- compute_lpd_features(series, config)
  lpd_features <- clean_features(lpd_features)
  fit <- fit_febama(lpd_features, config)
  weights <- prepare_febama_weights(fit, config, burnin = 0)
  forecast <- forecast_febama(series, config, lpd_features, fit = fit)
  performance <- summarize_performance(forecast)

  expect_equal(ncol(lpd_features$lpd), length(config$fore_model))
  expect_true(all(config$features_used[[1]] %in% colnames(lpd_features$feat)))
  expect_length(fit$beta, 1)
  expect_equal(dim(fit$beta[[1]]), c(1, 3))
  expect_length(weights, 1)
  expect_equal(length(weights[[1]]$beta), 3)
  expect_equal(dim(forecast$ff_feature), c(1, config$forecast_h))
  expect_equal(dim(forecast$w_time_varying), c(length(config$fore_model), config$forecast_h))
  expect_equal(dim(forecast$err_feature), c(1, 3))
  expect_equal(dim(performance), c(1, 3))
  expect_true(all(is.finite(forecast$ff_feature)))
  expect_true(all(is.finite(forecast$w_time_varying)))
})

test_that("compute_lpd_features supports the base parallel backend", {
  skip_if(parallel::detectCores() < 2, "parallel backend requires at least two cores")

  values <- tiny_series(36)
  config <- tiny_config(lpd_features_parallel = TRUE, ncores = 2)
  lpd_features <- compute_lpd_features(list(x = values), config)

  expect_equal(ncol(lpd_features$lpd), length(config$fore_model))
  expect_equal(ncol(lpd_features$feat), 42)
  expect_true(all(is.finite(lpd_features$lpd)))
})
