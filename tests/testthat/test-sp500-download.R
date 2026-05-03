source_sp500_downloader <- function() {
  script <- system.file("examples", "download_sp500.R", package = "febama")
  if (!nzchar(script)) {
    script <- testthat::test_path("..", "..", "inst", "examples", "download_sp500.R")
  }
  source(script, local = parent.frame())
}

test_that("S&P 500 downloader transformation computes paper-style returns", {
  source_sp500_downloader()

  prices <- data.frame(
    date = as.Date(c("2009-12-31", "2010-01-04", "2010-01-05")),
    close = c(100, 105, 103)
  )

  returns <- to_percent_log_returns(prices)

  expect_equal(returns$date, as.Date(c("2010-01-04", "2010-01-05")))
  expect_equal(returns$close, c(105, 103))
  expect_equal(
    returns$log_return,
    100 * diff(log(prices$close)),
    tolerance = 1e-12
  )
})

test_that("Yahoo URL uses the paper date window", {
  source_sp500_downloader()

  url <- yahoo_chart_url()

  expect_match(url, "query1[.]finance[.]yahoo[.]com", fixed = FALSE)
  expect_match(url, "%5EGSPC", fixed = TRUE)
  expect_match(url, paste0("period1=", date_to_unix(as.Date("2009-12-31"))), fixed = TRUE)
  expect_match(url, paste0("period2=", date_to_unix(as.Date("2019-09-19"))), fixed = TRUE)
  expect_match(url, "interval=1d", fixed = TRUE)
})
