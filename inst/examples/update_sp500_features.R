# Update rolling S&P 500 feature files for the paper example.
#
# Usage:
#   Rscript inst/examples/update_sp500_features.R data/sp500_daily_percent_log_returns.csv
#
# Default outputs:
#   data/sp500_features_all.csv
#   data/sp500_features_table3.csv
#
# The defaults match Section 4 of docs/paper.pdf:
#   - daily percent log returns;
#   - first forecasting origin after 1,250 returns;
#   - one-step forecast target;
#   - 100-day sliding window for feature calculation;
#   - Table 3's 15 selected stock-market features.

source_sp500_example <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  script <- if (length(file_arg) > 0L) {
    normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE)
  } else {
    system.file("examples", "update_sp500_features.R", package = "febama")
  }

  if (!nzchar(script)) {
    script <- file.path("inst", "examples", "update_sp500_features.R")
  }

  source(file.path(dirname(script), "sp500.R"), local = parent.frame())
}

source_sp500_example()

read_sp500_return_frame <- function(path) {
  if (!file.exists(path)) {
    stop("Cannot find S&P 500 return file: ", path, call. = FALSE)
  }

  data <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  normalized_names <- normalize_column_names(names(data))
  date_col <- first_matching_column(normalized_names, c("date", "time", "timestamp"))
  return_col <- first_matching_column(
    normalized_names,
    c("return", "returns", "log_return", "log_returns", "ret", "sp500_return")
  )

  if (is.na(return_col)) {
    stop("Input CSV must contain a daily percent log-return column.", call. = FALSE)
  }

  dates <- if (is.na(date_col)) {
    seq_along(data[[return_col]])
  } else {
    as.Date(data[[date_col]])
  }
  returns <- as.numeric(data[[return_col]])
  keep <- is.finite(returns) & !is.na(dates)

  out <- data.frame(
    date = dates[keep],
    log_return = returns[keep]
  )
  out <- out[order(out$date), , drop = FALSE]
  rownames(out) <- NULL
  out
}

compute_sp500_feature_frame <- function(returns,
                                        dates = seq_along(returns),
                                        history_burn = 1250L,
                                        forecast_h = 1L,
                                        feature_window = 100L,
                                        frequency = 1L) {
  if (length(returns) <= history_burn + forecast_h) {
    stop("Not enough returns for the requested history and forecast windows.",
         call. = FALSE)
  }

  tha_features_fun <- get("tha_features", envir = asNamespace("febama"))
  origins <- seq.int(history_burn, length(returns) - forecast_h)

  feature_rows <- lapply(origins, function(origin) {
    window_start <- max(1L, origin - feature_window + 1L)
    feature_input <- list(list(
      x = stats::ts(returns[seq.int(window_start, origin)], frequency = frequency)
    ))
    data.matrix(tha_features_fun(feature_input)[[1L]]$features)
  })
  feature_matrix <- do.call(rbind, feature_rows)

  data.frame(
    origin_index = origins,
    origin_date = dates[origins],
    forecast_index = origins + forecast_h,
    forecast_date = dates[origins + forecast_h],
    actual_log_return = returns[origins + forecast_h],
    feature_matrix,
    check.names = FALSE
  )
}

write_csv <- function(data, output) {
  output_dir <- dirname(output)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  utils::write.csv(data, output, row.names = FALSE)
  normalizePath(output, mustWork = FALSE)
}

write_feature_files <- function(features,
                                all_output,
                                selected_output,
                                selected_features = sp500_feature_set()) {
  missing_features <- setdiff(selected_features, colnames(features))
  if (length(missing_features) > 0L) {
    stop(
      "Missing selected S&P 500 features: ",
      paste(missing_features, collapse = ", "),
      call. = FALSE
    )
  }

  write_csv(features, all_output)

  id_cols <- c(
    "origin_index", "origin_date", "forecast_index",
    "forecast_date", "actual_log_return"
  )
  selected <- features[c(id_cols, selected_features)]
  write_csv(selected, selected_output)

  invisible(list(all = all_output, selected = selected_output))
}

parse_args <- function(args) {
  options <- list(
    input = if (length(args) > 0L && !startsWith(args[[1L]], "--")) {
      args[[1L]]
    } else {
      file.path("data", "sp500_daily_percent_log_returns.csv")
    },
    all_output = file.path("data", "sp500_features_all.csv"),
    selected_output = file.path("data", "sp500_features_table3.csv"),
    history_burn = 1250L,
    forecast_h = 1L,
    feature_window = 100L
  )

  named_args <- args[startsWith(args, "--")]
  for (arg in named_args) {
    parts <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1L]]
    if (length(parts) != 2L) {
      stop("Arguments must use --name=value syntax: ", arg, call. = FALSE)
    }

    name <- gsub("-", "_", parts[[1L]])
    value <- parts[[2L]]
    if (!name %in% names(options)) {
      stop("Unknown argument: --", parts[[1L]], call. = FALSE)
    }

    options[[name]] <- switch(
      name,
      history_burn = as.integer(value),
      forecast_h = as.integer(value),
      feature_window = as.integer(value),
      value
    )
  }

  options
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  options <- parse_args(args)
  data <- read_sp500_return_frame(options$input)
  features <- compute_sp500_feature_frame(
    returns = data$log_return,
    dates = data$date,
    history_burn = options$history_burn,
    forecast_h = options$forecast_h,
    feature_window = options$feature_window
  )
  write_feature_files(
    features = features,
    all_output = options$all_output,
    selected_output = options$selected_output
  )

  message("Wrote ", nrow(features), " rolling S&P 500 feature rows.")
  message("All features: ", normalizePath(options$all_output, mustWork = FALSE))
  message("Table 3 features: ", normalizePath(options$selected_output, mustWork = FALSE))
  invisible(features)
}

if (sys.nframe() == 0L) {
  main()
}
