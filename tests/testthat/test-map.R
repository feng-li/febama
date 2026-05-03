map_training_data <- function(n_obs = 8, n_models = 2, n_features = 3) {
  lpd <- matrix(log(seq(0.2, 0.8, length.out = n_obs * n_models)),
                nrow = n_obs, ncol = n_models)
  feat <- matrix(seq(-1, 1, length.out = n_obs * n_features),
                 nrow = n_obs, ncol = n_features)
  colnames(lpd) <- paste0("m", seq_len(n_models))
  colnames(feat) <- paste0("f", seq_len(n_features))
  list(lpd = lpd, feat = feat)
}

map_config <- function(n_models = 2,
                       variable_selection = FALSE,
                       varsel_init = "all-in",
                       n_iter = 2,
                       init_optim = FALSE) {
  febama_config(
    fore_model = paste0("m", seq_len(n_models)),
    features_used = c("f1", "f2", "f3"),
    variable_selection = variable_selection,
    varsel_init = varsel_init,
    alg_name = "MAP",
    n_iter = n_iter,
    init_optim = init_optim,
    beta_shrinkage = 100
  )
}

test_that("MAP keeps all features active when variable selection is disabled", {
  set.seed(1)
  config <- map_config(variable_selection = FALSE, varsel_init = "all-out")
  data <- map_training_data()

  fit <- fit_febama(data, config)

  expect_equal(config$varSelArgs[[1]]$init, "all-in")
  expect_true(all(fit$betaIdx[[1]] == 1))
  expect_true(all(is.finite(fit$beta[[1]])))
})

test_that("MAP respects explicit variable-selection candidate subsets", {
  set.seed(2)
  config <- map_config(variable_selection = TRUE, varsel_init = "all-out", n_iter = 3)
  config$varSelArgs[[1]]$cand <- 3L
  data <- map_training_data()

  fit <- fit_febama(data, config)

  expect_equal(febama:::resolve_varsel_candidates(4, c(1L, 3L, 3L)), 3L)
  expect_equal(febama:::initialize_beta_idx(4, config$varSelArgs[[1]]), c(1, 1, 0, 1))
  expect_true(all(fit$betaIdx[[1]][, c(1, 2, 4)] == 1))
})

test_that("MAP optimizer only updates active coefficients", {
  data <- map_training_data()
  config <- map_config(n_iter = 1)
  beta <- list(c(0.1, 99, -0.2, 0.3))
  beta_idx <- list(c(1, 0, 1, 1))

  out <- febama:::optimize_map_coefficients(
    data = data,
    beta = beta,
    betaIdx = beta_idx,
    priArgs = config$priArgs,
    varSelArgs = config$varSelArgs,
    features_used = config$features_used,
    model_update = 1
  )

  expect_equal(out$beta[[1]][2], 0)
  expect_true(all(is.finite(out$beta[[1]][c(1, 3, 4)])))
})

test_that("logscore remains finite for very small predictive densities", {
  data <- list(
    lpd = matrix(c(-1000, -1001, -1002, -1e6, -1e6 - 1, -1e6 - 2),
                 nrow = 2, byrow = TRUE),
    feat = matrix(c(0, 1), ncol = 1)
  )
  colnames(data$feat) <- "f1"
  beta <- list(c(0, 0), c(0, 0))
  beta_idx <- list(c(1, 1), c(1, 1))
  features_used <- list("f1", "f1")

  scores <- logscore(data, beta, beta_idx, features_used, sum = FALSE)

  expect_true(all(is.finite(scores)))
  expect_lt(scores[[2]], -999999)
})

test_that("Bernoulli betaIdx prior reads component prior probability", {
  pri <- febama:::default_prior_args(beta_idx_type = "bern")
  pri$betaIdx$prob <- 0.25

  val <- log_priors(
    beta = list(c(0, 0, 0)),
    betaIdx = list(c(1, 1, 0)),
    varSelArgs = list(list(cand = 2:3, init = "all-in")),
    priArgs = list(pri),
    sum = TRUE
  )

  expect_true(is.finite(val))
})
