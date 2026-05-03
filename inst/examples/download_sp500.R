# Download and transform S&P 500 data for the paper example.
#
# Usage:
#   Rscript inst/examples/download_sp500.R data/sp500_daily_percent_log_returns.csv
#
# The paper uses daily percent log returns from January 4, 2010 to
# September 18, 2019. This script downloads S&P 500 daily closes from Yahoo
# Finance's chart endpoint starting on the previous trading day, then computes
#
#   log_return = 100 * diff(log(close))
#
# Output columns:
#   date, close, log_return
#
# The output CSV is accepted by inst/examples/sp500.R.

sp500_paper_start <- as.Date("2010-01-04")
sp500_paper_end <- as.Date("2019-09-18")
sp500_download_start <- as.Date("2009-12-31")

date_to_unix <- function(date) {
  as.integer(as.POSIXct(as.Date(date), tz = "UTC"))
}

yahoo_chart_url <- function(symbol = "^GSPC",
                            from = sp500_download_start,
                            to = sp500_paper_end) {
  query <- paste0(
    "period1=", date_to_unix(from),
    "&period2=", date_to_unix(as.Date(to) + 1L),
    "&interval=1d",
    "&events=history",
    "&includeAdjustedClose=true"
  )
  paste0(
    "https://query1.finance.yahoo.com/v8/finance/chart/",
    utils::URLencode(symbol, reserved = TRUE),
    "?",
    query
  )
}

download_yahoo_daily <- function(symbol = "^GSPC",
                                 from = sp500_download_start,
                                 to = sp500_paper_end,
                                 destfile = tempfile(fileext = ".json")) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("The downloader requires the 'jsonlite' package.", call. = FALSE)
  }

  url <- yahoo_chart_url(symbol = symbol, from = from, to = to)
  status <- utils::download.file(
    url = url,
    destfile = destfile,
    mode = "wb",
    quiet = TRUE
  )
  if (!identical(status, 0L)) {
    stop("Failed to download S&P 500 data from Yahoo Finance.", call. = FALSE)
  }

  payload <- jsonlite::fromJSON(destfile, simplifyVector = FALSE)
  error <- payload$chart$error
  if (!is.null(error)) {
    stop("Yahoo Finance returned an error: ", error$description, call. = FALSE)
  }

  result <- payload$chart$result[[1L]]
  timestamps <- unlist(result$timestamp, use.names = FALSE)
  quote <- result$indicators$quote[[1L]]
  close <- unlist(quote$close, use.names = FALSE)
  adjclose <- result$indicators$adjclose[[1L]]$adjclose
  if (!is.null(adjclose)) {
    adjclose <- unlist(adjclose, use.names = FALSE)
    if (length(adjclose) == length(close) && any(is.finite(adjclose))) {
      close <- adjclose
    }
  }

  prices <- data.frame(
    date = as.Date(as.POSIXct(timestamps, origin = "1970-01-01", tz = "UTC")),
    close = as.numeric(close)
  )
  prices <- prices[order(prices$date), , drop = FALSE]
  prices <- prices[is.finite(prices$close) & !is.na(prices$date), , drop = FALSE]
  rownames(prices) <- NULL

  if (nrow(prices) < 2L) {
    stop("Downloaded price data has fewer than two observations.", call. = FALSE)
  }

  prices
}

to_percent_log_returns <- function(prices,
                                   start = sp500_paper_start,
                                   end = sp500_paper_end) {
  prices <- prices[order(prices$date), , drop = FALSE]
  transformed <- data.frame(
    date = prices$date[-1L],
    close = prices$close[-1L],
    log_return = 100 * diff(log(prices$close))
  )
  transformed <- transformed[
    transformed$date >= as.Date(start) & transformed$date <= as.Date(end),
    ,
    drop = FALSE
  ]
  rownames(transformed) <- NULL

  if (nrow(transformed) == 0L) {
    stop("No transformed returns are available in the requested paper window.",
         call. = FALSE)
  }
  if (!identical(transformed$date[[1L]], as.Date(start))) {
    warning(
      "First transformed return is ", transformed$date[[1L]],
      ", not the paper start date ", as.Date(start),
      ". Check whether the source includes the prior trading day.",
      call. = FALSE
    )
  }

  transformed
}

write_csv <- function(data, output) {
  output_dir <- dirname(output)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  utils::write.csv(data, output, row.names = FALSE)
  normalizePath(output, mustWork = FALSE)
}

parse_args <- function(args) {
  options <- list(
    output = if (length(args) > 0L && !startsWith(args[[1L]], "--")) {
      args[[1L]]
    } else {
      file.path("data", "sp500_daily_percent_log_returns.csv")
    },
    symbol = "^GSPC",
    start = sp500_paper_start,
    end = sp500_paper_end,
    download_start = sp500_download_start,
    raw_output = NA_character_
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
      start = as.Date(value),
      end = as.Date(value),
      download_start = as.Date(value),
      value
    )
  }

  options
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  options <- parse_args(args)
  prices <- download_yahoo_daily(
    symbol = options$symbol,
    from = options$download_start,
    to = options$end
  )
  returns <- to_percent_log_returns(
    prices,
    start = options$start,
    end = options$end
  )

  if (!is.na(options$raw_output)) {
    write_csv(prices, options$raw_output)
  }

  output <- write_csv(returns, options$output)
  message("Wrote ", nrow(returns), " S&P 500 daily percent log returns to ", output)
  message("Date range: ", min(returns$date), " to ", max(returns$date))
  invisible(returns)
}

if (sys.nframe() == 0L) {
  main()
}
