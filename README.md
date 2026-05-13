# febama

[![R-CMD-check](https://github.com/feng-li/febama/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/feng-li/febama/actions/workflows/R-CMD-check.yaml)


`febama` is a framework for density forecast combination with time-varying
weights based on time-series features. Each base model contributes a predictive
density, features map to model weights through a softmax-style regression, and
the coefficient vectors are estimated with Bayesian log predictive scores and
optional variable selection.

## Python Port

The FEBAMA workflow has been ported to Python as [`gsm.febama`](https://github.com/feng-li/gsm), a submodule of the
general smooth-mixture package `gsm`. The Python port keeps the FEBAMA idea of feature-driven Bayesian forecast
averaging, but places it inside the broader GSM/MoE codebase with JAX-based scoring, standard Python data tooling, and
pluggable predictive distributions.

New Python development should target the `gsm.febama` module. This R package
remains the original implementation and a reference for API behavior, S&P 500
examples, feature construction, and paper replication checks.

In the Python `gsm` repository, install and run the current FEBAMA example
with:

```sh
python -m pip install -e ".[dev,febama]"
python scripts/run_febama_example.py --max-origins 4 --test-size 1 --max-iter 100
```

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

The core package uses `forecast`, `mvtnorm`, `parallel`, and `tsfeatures`.
The paper-style S&P 500 volatility example additionally needs optional
packages:

```r
install.packages(c("rugarch", "highfrequency", "xts", "stochvol"))
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

## S&P 500 Paper Example

The paper's stock-market experiment is implemented as reproducible scripts in
`inst/examples/`. By default these scripts use the paper-style settings:
one-step forecasts, a 1,250-trading-day rolling model window, a 100-day feature
window, GARCH/RGARCH/SV base models, and the 15 stock-market features from
Table 3.

The generated S&P 500 CSV files are written under `data/` and are ignored by
git. Recreate them locally with the following steps.

Download and transform the daily S&P 500 closes from Yahoo Finance into the
paper's daily percent log returns:

```sh
Rscript inst/examples/download_sp500.R data/sp500_daily_percent_log_returns.csv
```

This writes `date`, `close`, and `log_return`, where `log_return` is
`100 * diff(log(close))` over January 4, 2010 to September 18, 2019.

Update the rolling feature files used by the S&P 500 experiment:

```sh
Rscript inst/examples/update_sp500_features.R data/sp500_daily_percent_log_returns.csv
```

This writes `data/sp500_features_all.csv` with all 42 THA features and
`data/sp500_features_table3.csv` with the 15 stock-market features listed in
Table 3 of the paper.

Run the paper-style FEBAMA example:

```sh
Rscript inst/examples/sp500.R data/sp500_daily_percent_log_returns.csv
```

Compare MAP and SGLD on the same rolling-origin S&P 500 density/features:

```sh
Rscript inst/examples/sp500_compare_algorithms.R data/sp500_daily_percent_log_returns.csv
```

Write the per-origin comparison table to CSV by passing an output path:

```sh
Rscript inst/examples/sp500_compare_algorithms.R \
  data/sp500_daily_percent_log_returns.csv \
  data/sp500_map_sgld_comparison.csv
```

For a plumbing-only run without the optional volatility-model packages:

```sh
FEBAMA_SP500_FAST=1 Rscript inst/examples/sp500.R data/sp500_daily_percent_log_returns.csv
FEBAMA_SP500_FAST=1 Rscript inst/examples/sp500_compare_algorithms.R data/sp500_daily_percent_log_returns.csv
```

Control the number of rolling origins with `FEBAMA_SP500_ORIGINS`. For example:

```sh
FEBAMA_SP500_ORIGINS=3 Rscript inst/examples/sp500.R data/sp500_daily_percent_log_returns.csv
```

The comparison script defaults to `n_iter = 2` so both algorithms perform an
outer inference update. SGLD settings can be adjusted with
`FEBAMA_SP500_SGLD_N_EPOCH`, `FEBAMA_SP500_SGLD_MAX_BATCH_SIZE`,
`FEBAMA_SP500_SGLD_BURNIN_PROP`, and `FEBAMA_SP500_SGLD_STEPSIZE`.

A one-origin fast-mode run on the bundled local S&P 500 return file produced:

| Algorithm | Origin | Forecast | Actual | Log score | MASE | SMAPE | Time |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| MAP | 2443 | 0.257851 | 0.034263 | -1.446465 | 0.233567 | 153.0823 | 0.398s |
| SGLD | 2443 | 0.257843 | 0.034263 | -1.446450 | 0.233559 | 153.0811 | 0.492s |

In that smoke comparison, SGLD was slightly better on log score, MASE, and
SMAPE, while MAP was slightly faster. The differences are tiny for this single
origin; use more origins and the paper-style volatility models for a substantive
empirical comparison.

## Model Summary

For each historical cutoff, FEBAMA:

1. fits each base forecaster in `config$fore_model`;
2. evaluates each model's Gaussian predictive density for the next holdout
   point;
3. computes the THA time-series feature set using `tsfeatures`;
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

inst/examples/          # S&P 500 download, feature, and forecast scripts
man/                    # generated package documentation
docs/                   # methodology notes
data/                   # saved training artifacts
tests/testthat/         # focused public API and S&P 500 helper tests
.github/workflows/      # R CMD check workflow
```

## Development

Run the local checks before pushing:

```sh
Rscript -e "pkgload::load_all('.', quiet = TRUE); testthat::test_dir('tests/testthat')"
R CMD build .
R CMD check febama_0.0.0.9000.tar.gz
```

The GitHub Actions workflow runs `R CMD check` on Linux, macOS, Windows, and
R-devel. Local S&P 500 data outputs are intentionally ignored so the repository
contains reproducible scripts rather than downloaded market data.
