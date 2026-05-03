# febama

Feature-based Bayesian Forecasting Model Averaging.

`febama` is a framework for density forecast combination with time-varying
weights based on time-series features. Each base model contributes a predictive
density, features map to model weights through a softmax-style regression, and
the coefficient vectors are estimated with Bayesian log predictive scores and
optional variable selection.

The method is based on the published paper

```
@article{LiL2023BayesianForecast,
  title = {Bayesian Forecast Combination Using Time-Varying Features},
  author = {Li, Li and Kang, Yanfei and Li, Feng},
  date = {2023-07},
  journaltitle = {International Journal of Forecasting},
  volume = {39},
  number = {3},
  pages = {1287--1302},
  issn = {0169-2070},
  doi = {10.1016/j.ijforecast.2022.06.002},
  url = {https://arxiv.org/abs/2108.02082},
  urldate = {2023-06-21},
  language = {en},
  keywords = {Bayesian density forecasting,Forecast combination,Interpretability,Log predictive score,Time-varying features}
}

```


## Installation

Install the package from GitHub with:

```r
devtools::install_github("feng-li/febama")
```

## Public API

The recommended package entry points are:

| Function | Purpose |
| --- | --- |
| `febama_config()` | Build the nested model configuration. |
| `compute_lpd_features()` | Build historical log predictive densities and feature matrices. |
| `clean_features()` | Drop unusable feature columns and keep scaling metadata. |
| `fit_febama()` | Estimate feature-weight coefficients. |
| `forecast_febama()` | Produce recursive combined forecasts. |
| `summarize_performance()` | Aggregate log score, MASE, and SMAPE. |

The lower-level functions, such as `lpd_features_multi()`, `febama_mcmc()`,
and `forecast_feature_results_multi()`, remain exported for compatibility with
the original research workflow.

## Data Shape

The training and forecasting functions expect each series to be a list with:

- `x`: the in-sample time series.
- `xx`: the held-out future values used by the current forecasting routine for
  scoring.

This matches the shape used by the M3/M4 examples in the original experiments.

## Basic Workflow

```r
library(febama)

config <- febama_config(
  frequency = 12,
  forecast_h = 18,
  train_h = 1,
  history_burn = 25,
  fore_model = c("ets_fore", "naive_fore", "rw_drift_fore", "auto.arima_fore"),
  features_used = c("x_acf1", "diff1_acf1", "entropy", "alpha", "beta", "unitroot_kpss")
)

lpd_features <- compute_lpd_features(series, config)
lpd_features <- clean_features(lpd_features)

fit <- fit_febama(lpd_features, config)

forecast <- forecast_febama(
  data = series,
  config = config,
  lpd_features = lpd_features,
  fit = fit
)

summarize_performance(forecast)
```

For a list of series, pass the list to the same functions:

```r
lpd_features <- compute_lpd_features(series_list, config)
lpd_features <- clean_features(lpd_features)
fits <- fit_febama(lpd_features, config)
forecasts <- forecast_febama(series_list, config, lpd_features, fit = fits)
summarize_performance(forecasts)
```

## Model Summary

For each historical cutoff, FEBAMA:

1. fits each base forecaster in `config$fore_model`;
2. evaluates each model's Gaussian predictive density for the next holdout
   point;
3. computes time-series features using `M4metalearning::THA_features()`;
4. estimates feature coefficients by maximizing or sampling the log predictive
   score with coefficient priors and optional feature-selection indicators.

With `m` base forecasters, FEBAMA estimates `m - 1` feature-weight regressions.
The last forecaster is the baseline component, which keeps the weights
identified and makes the full weight vector sum to one.

## Repository Layout

```text
R/
  api.R                 # public wrappers and configuration helpers
  default_parameters.R  # default model settings
  densities_features.R  # historical density and feature construction
  models.R              # base forecasting model wrappers
  logscore.R            # log predictive score and gradient
  priors.R              # coefficient and feature-selection priors
  posterior.R           # posterior objective and gradient
  mcmc.R                # MAP/SGLD/MCMC fitting
  forecast.R            # recursive forecasting and performance summaries

man/                    # generated package documentation
docs/                   # methodology notes
data/                   # saved training artifacts
```
