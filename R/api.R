#' Create a FEBAMA model configuration
#'
#' This constructor keeps the package examples and scripts from hand-building
#' large nested configuration lists.
#'
#' @param frequency Seasonal frequency of the input series.
#' @param ets_model ETS model code kept for compatibility with older scripts.
#' @param ts_scale Should the time series be standardized before fitting base forecasters?
#' @param forecast_h Forecast horizon.
#' @param train_h Historical holdout horizon used when building training scores.
#' @param history_burn Minimum history length before rolling training starts.
#' @param PI_level Prediction interval level used to infer forecast standard deviations.
#' @param roll Optional rolling-history window for fitting base forecasters.
#' @param feature_window Optional rolling-history window for feature calculation.
#' @param fore_model Character vector of base forecaster function names.
#' @param features_used Character vector shared by each non-baseline model, or a
#'   list with one character vector per non-baseline model.
#' @param lpd_features_parallel Should historical predictive densities and
#'   features be computed in parallel?
#' @param ncores Number of worker cores to use when `lpd_features_parallel` is
#'   `TRUE`.
#' @param variable_selection Should feature-selection indicators be sampled?
#' @param varsel_init Initial feature-selection state.
#' @param beta_shrinkage Diagonal covariance scale for the coefficient prior.
#' @param init_optim Should BFGS be used to initialize coefficients?
#' @param alg_name Inference algorithm, either `"MAP"` or `"sgld"`.
#' @param n_iter Number of outer MCMC iterations.
#' @param max_batch_size Maximum SGLD minibatch size.
#' @param n_epoch Number of SGLD epochs inside each outer update.
#' @param burnin_prop Within-SGLD burn-in proportion for Polyak averaging.
#' @param stepsize Fixed SGLD step size. Use `NA` for the decaying schedule.
#' @param sgld_gamma Exponent for the decaying SGLD step-size schedule.
#' @param sgld_a Scale for the decaying SGLD step-size schedule.
#' @param sgld_b Offset for the decaying SGLD step-size schedule.
#'
#' @return A list consumable by the FEBAMA training and forecasting functions.
#' @export
febama_config <- function(frequency = 12,
                          ets_model = "AAN",
                          ts_scale = TRUE,
                          forecast_h = 18,
                          train_h = 1,
                          history_burn = 25,
                          PI_level = 90,
                          roll = NULL,
                          feature_window = NULL,
                          fore_model = c("ets_fore", "naive_fore", "rw_drift_fore", "auto.arima_fore"),
                          features_used = c("x_acf1", "diff1_acf1", "entropy", "alpha", "beta", "unitroot_kpss"),
                          lpd_features_parallel = FALSE,
                          ncores = 1,
                          variable_selection = TRUE,
                          varsel_init = c("all-in", "all-out", "random"),
                          beta_shrinkage = 10,
                          init_optim = TRUE,
                          alg_name = c("MAP", "sgld"),
                          n_iter = 1,
                          max_batch_size = 108,
                          n_epoch = 10,
                          burnin_prop = 0.4,
                          stepsize = 0.1,
                          sgld_gamma = 0.55,
                          sgld_a = 0.4,
                          sgld_b = 10) {
    varsel_init <- match.arg(varsel_init)
    alg_name <- match.arg(alg_name)

    if (length(fore_model) < 2) {
        stop("`fore_model` must contain at least two forecasting functions.")
    }

    num_models_updated <- length(fore_model) - 1
    features_by_model <- normalize_features_used(features_used, num_models_updated)
    varsel_cand <- if (isTRUE(variable_selection)) "2:end" else NULL

    list(
        frequency = frequency,
        ets_model = ets_model,
        ts_scale = ts_scale,
        forecast_h = forecast_h,
        train_h = train_h,
        history_burn = history_burn,
        PI_level = PI_level,
        roll = roll,
        feature_window = feature_window,
        features_used = features_by_model,
        fore_model = fore_model,
        lpd_features_parl = list(par = lpd_features_parallel, ncores = ncores),
        varSelArgs = rep(list(list(cand = varsel_cand, init = varsel_init)), num_models_updated),
        priArgs = rep(list(default_prior_args(shrinkage = beta_shrinkage)), num_models_updated),
        algArgs = default_alg_args(
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
    )
}

#' Compute historical predictive densities and features
#'
#' @param data A single series list with element `x`, or a list of such series.
#' @param config A FEBAMA configuration created by [febama_config()].
#'
#' @return A training-data list with `lpd` and `feat`, or a list of those objects.
#' @export
compute_lpd_features <- function(data, config) {
    if (is_series_data(data)) {
        return(lpd_features_multi(data, config))
    }

    lapply(data, lpd_features_multi, model_conf = config)
}

#' Clean feature matrices
#'
#' Removes feature columns containing missing values and stores the scaling
#' parameters needed for forecasting.
#'
#' @param lpd_features A single training-data object or a list of them.
#'
#' @return Cleaned training data in the same single/list shape as the input.
#' @export
clean_features <- function(lpd_features) {
    if (is_lpd_feature_data(lpd_features)) {
        return(feature_clean(list(lpd_features))[[1]])
    }

    feature_clean(lpd_features)
}

#' Fit a FEBAMA model
#'
#' @param lpd_features A single cleaned training-data object or a list of them.
#' @param config A FEBAMA configuration created by [febama_config()].
#'
#' @return MCMC output from [febama_mcmc()], or a list of such outputs.
#' @export
fit_febama <- function(lpd_features, config) {
    if (is_lpd_feature_data(lpd_features)) {
        return(febama_mcmc(lpd_features, config))
    }

    lapply(lpd_features, febama_mcmc, model_conf = config)
}

#' Prepare posterior coefficient summaries for forecasting
#'
#' @param fit MCMC output from [fit_febama()], or a list of fit objects.
#' @param config A FEBAMA configuration created by [febama_config()].
#' @param burnin Burn-in as a proportion in `[0, 1)` or an iteration count.
#' @param inclusion_threshold Posterior inclusion probability required for a
#'   feature to be used in the forecast weights.
#' @param statistic Posterior summary statistic for coefficients.
#'
#' @return A `beta_pre` object consumed by the legacy forecasting routine, or a
#'   list of such objects.
#' @export
prepare_febama_weights <- function(fit,
                                   config,
                                   burnin = 0.5,
                                   inclusion_threshold = 0.5,
                                   statistic = c("mean", "median")) {
    statistic <- match.arg(statistic)

    if (!is_febama_fit(fit)) {
        return(lapply(
            fit,
            prepare_febama_weights,
            config = config,
            burnin = burnin,
            inclusion_threshold = inclusion_threshold,
            statistic = statistic
        ))
    }

    stat_fun <- switch(
        statistic,
        mean = colMeans,
        median = function(x) apply(x, 2, stats::median)
    )

    out <- vector("list", length(fit$beta))
    for (iComp in seq_along(fit$beta)) {
        beta_samples <- fit$beta[[iComp]]
        beta_idx_samples <- fit$betaIdx[[iComp]]
        rows <- posterior_rows(nrow(beta_samples), burnin)

        inclusion_prob <- colMeans(beta_idx_samples[rows, , drop = FALSE])
        selected <- inclusion_prob >= inclusion_threshold
        selected[1] <- TRUE

        selected_feature_pos <- which(selected[-1])
        feature_names <- config$features_used[[iComp]]
        features_select <- feature_names[selected_feature_pos]
        if (is.null(features_select)) {
            features_select <- selected_feature_pos
        }

        beta <- stat_fun(beta_samples[rows, selected, drop = FALSE])
        names(beta) <- c("(Intercept)", features_select)

        out[[iComp]] <- list(
            beta = as.numeric(beta),
            features_select = features_select,
            inclusion_prob = inclusion_prob
        )
    }

    out
}

#' Forecast with a fitted FEBAMA model
#'
#' @param data A single series list with elements `x` and `xx`, or a list of such
#'   series. The current implementation uses `xx` to score the forecast.
#' @param config A FEBAMA configuration created by [febama_config()].
#' @param lpd_features Cleaned training data matching `data`.
#' @param fit MCMC output from [fit_febama()].
#' @param intercept Should an intercept be included in the weight model?
#'
#' @return A forecast result, or a list of forecast results.
#' @export
forecast_febama <- function(data,
                            config,
                            lpd_features,
                            fit,
                            intercept = TRUE) {
    if (missing(fit) || is.null(fit)) {
        stop("`fit` must be provided.")
    }

    if (is_series_data(data)) {
        return(forecast_feature_results_multi(
            ts = data,
            model_conf = config,
            intercept = intercept,
            data = lpd_features,
            beta_out = fit
        ))
    }

    mapply(
        function(series, features, model_fit) {
            forecast_feature_results_multi(
                ts = series,
                model_conf = config,
                intercept = intercept,
                data = features,
                beta_out = model_fit
            )
        },
        series = data,
        features = lpd_features,
        model_fit = fit,
        SIMPLIFY = FALSE,
        USE.NAMES = FALSE
    )
}

#' Summarize FEBAMA forecast performance
#'
#' @param forecasts A forecast result from [forecast_febama()] or a list of them.
#'
#' @return Matrix with log score, MASE, and SMAPE.
#' @export
summarize_performance <- function(forecasts) {
    if (is_series_data(forecasts) && !is.null(forecasts$err_feature)) {
        forecasts <- list(forecasts)
    }

    forecast_feature_performance(forecasts)
}

default_prior_args <- function(mean = 0,
                               covariance = "identity",
                               shrinkage = 1,
                               beta_idx_type = "beta",
                               alpha0 = 1,
                               beta0 = 1) {
    list(
        beta = list(
            type = "cond-mvnorm",
            mean = mean,
            covariance = covariance,
            shrinkage = shrinkage
        ),
        betaIdx = list(
            type = beta_idx_type,
            alpha0 = alpha0,
            beta0 = beta0
        )
    )
}

default_alg_args <- function(init_optim,
                             alg_name,
                             n_iter,
                             max_batch_size,
                             n_epoch,
                             burnin_prop,
                             stepsize,
                             sgld_gamma,
                             sgld_a,
                             sgld_b) {
    list(
        initOptim = init_optim,
        algName = alg_name,
        nIter = n_iter,
        sgld = list(
            max_batchSize = max_batch_size,
            nEpoch = n_epoch,
            burninProp = burnin_prop,
            stepsize = stepsize,
            gama = sgld_gamma,
            a = sgld_a,
            b = sgld_b
        )
    )
}

normalize_features_used <- function(features_used, num_models_updated) {
    if (is.null(features_used) || length(features_used) == 0) {
        return(rep(list(NULL), num_models_updated))
    }

    if (is.list(features_used) && !is.data.frame(features_used)) {
        if (length(features_used) != num_models_updated) {
            stop("`features_used` must have one entry per non-baseline model.")
        }
        return(features_used)
    }

    rep(list(features_used), num_models_updated)
}

posterior_rows <- function(n_iter, burnin) {
    if (!is.numeric(burnin) || length(burnin) != 1 || is.na(burnin) || burnin < 0) {
        stop("`burnin` must be a non-negative number.")
    }

    n_drop <- if (burnin < 1) {
        floor(n_iter * burnin)
    } else {
        as.integer(burnin)
    }
    n_drop <- min(n_drop, n_iter - 1)

    seq.int(n_drop + 1, n_iter)
}

is_series_data <- function(x) {
    is.list(x) && !is.null(x$x)
}

is_lpd_feature_data <- function(x) {
    is.list(x) && !is.null(x$lpd) && !is.null(x$feat)
}

is_febama_fit <- function(x) {
    is.list(x) && !is.null(x$beta) && !is.null(x$betaIdx)
}
