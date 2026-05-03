sgld_training_data <- function(n_obs = 10, n_models = 3, n_features = 2) {
  lpd <- matrix(log(seq(0.2, 0.8, length.out = n_obs * n_models)),
                nrow = n_obs, ncol = n_models)
  feat <- matrix(seq(-1, 1, length.out = n_obs * n_features),
                 nrow = n_obs, ncol = n_features)
  colnames(lpd) <- paste0("m", seq_len(n_models))
  colnames(feat) <- paste0("f", seq_len(n_features))
  list(lpd = lpd, feat = feat)
}

sgld_config <- function(n_models = 3,
                        variable_selection = FALSE,
                        n_epoch = 2,
                        max_batch_size = 4,
                        burnin_prop = 0,
                        stepsize = 1e-4) {
  febama_config(
    fore_model = paste0("m", seq_len(n_models)),
    features_used = c("f1", "f2"),
    variable_selection = variable_selection,
    varsel_init = "all-in",
    alg_name = "sgld",
    n_iter = 2,
    init_optim = FALSE,
    max_batch_size = max_batch_size,
    n_epoch = n_epoch,
    burnin_prop = burnin_prop,
    stepsize = stepsize,
    beta_shrinkage = 100
  )
}

test_that("SGLD honors configured epochs and max batch size", {
  set.seed(1)
  data <- sgld_training_data(n_obs = 10)
  config <- sgld_config(n_epoch = 2, max_batch_size = 4)
  beta <- rep(list(c(0, 0, 0)), 2)
  beta_idx <- rep(list(c(1, 1, 1)), 2)

  out <- febama:::SGLD_gibbs(data, beta, beta_idx, config)

  expect_equal(nrow(out$beta_sgld[[1]]), 2 * ceiling(10 / 4))
  expect_equal(nrow(out$beta_sgld[[2]]), 2 * ceiling(10 / 4))
  expect_length(out$accept_prob, 2)
  expect_equal(unlist(out$accept_prob), c(1, 1))

  batches <- febama:::make_sgld_batches(10, 4)
  expect_equal(length(batches), 3)
  expect_lte(max(lengths(batches)), 4)
  expect_setequal(unlist(batches, use.names = FALSE), seq_len(10))
})

test_that("SGLD decaying step sizes keep decaying", {
  args <- list(a = 0.4, b = 10, gama = 0.55)

  decayed <- vapply(
    seq_len(4),
    function(i) febama:::sgld_step_size(NA_real_, i, args),
    numeric(1)
  )
  fixed <- vapply(
    seq_len(4),
    function(i) febama:::sgld_step_size(0.01, i, args),
    numeric(1)
  )

  expect_true(all(diff(decayed) < 0))
  expect_equal(fixed, rep(0.01, 4))
})

test_that("SGLD burn-in validation handles edge cases", {
  data <- sgld_training_data(n_obs = 10, n_models = 2)
  beta <- list(c(0, 0, 0))
  beta_idx <- list(c(1, 1, 1))

  set.seed(2)
  config <- sgld_config(n_models = 2, n_epoch = 1, max_batch_size = 20, burnin_prop = 0)
  out <- febama:::SGLD_gibbs(data, beta, beta_idx, config)
  expect_equal(nrow(out$beta_sgld[[1]]), 1)
  expect_true(all(is.finite(out$beta[[1]])))

  config$algArgs$sgld$burninProp <- 1
  expect_error(
    febama:::SGLD_gibbs(data, beta, beta_idx, config),
    "burninProp"
  )
})

test_that("public SGLD fit returns per-component diagnostics", {
  set.seed(3)
  data <- sgld_training_data(n_obs = 10, n_models = 3)
  config <- sgld_config(variable_selection = TRUE, n_epoch = 2, max_batch_size = 4)

  fit <- fit_febama(data, config)

  expect_length(fit$beta, 2)
  expect_length(fit$betaIdx, 2)
  expect_length(fit$accept_prob, 2)
  expect_equal(dim(fit$beta[[1]]), c(2, 3))
  expect_equal(dim(fit$beta[[2]]), c(2, 3))
  expect_equal(dim(fit$accept_prob[[1]]), c(2, 1))
  expect_equal(dim(fit$accept_prob[[2]]), c(2, 1))
  expect_true(all(is.finite(fit$accept_prob[[1]])))
  expect_true(all(is.finite(fit$accept_prob[[2]])))
})
