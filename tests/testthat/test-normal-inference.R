# Moment-only inference tests: no PageRank or permutation jobs are needed.

test_that("normal inference preserves ordinary Gaussian upper-tail probabilities", {
  beta = c(-8, -3, -1, 0, 1, 3, 8)
  null_mean = seq_along(beta) / 100
  null_sd = seq_along(beta) / 1000
  observed = null_mean + beta * null_sd
  result = tkoi_normal_inference(observed, null_mean, null_sd)

  expect_s3_class(result, "data.frame")
  expect_named(result, c("beta", "p_value", "log_p_value", "p_value_bounded", "inference_status"))
  expect_equal(result$beta, beta, tolerance = 1e-13)
  expect_equal(result$log_p_value, stats::pnorm(result$beta, lower.tail = FALSE, log.p = TRUE), tolerance = 0)
  expect_equal(result$p_value, stats::pnorm(result$beta, lower.tail = FALSE), tolerance = 1e-14)
  expect_identical(result$p_value_bounded, rep(FALSE, length(beta)))
  expect_identical(result$inference_status, rep("testable", length(beta)))
})

test_that("extreme finite beta retains Gaussian log probabilities and flags bounded numeric values", {
  beta = c(35, 38, 40, 100, 1000, 1e6)
  result = tkoi_normal_inference(beta, rep(0, length(beta)), rep(1, length(beta)))
  log_expected = stats::pnorm(beta, lower.tail = FALSE, log.p = TRUE)
  floor_log = log(.Machine$double.xmin)

  expect_identical(result$beta, beta)
  expect_equal(result$log_p_value, log_expected, tolerance = 0)
  expect_true(all(is.finite(result$log_p_value)))
  expect_true(all(diff(result$log_p_value) < 0))
  expect_identical(result$p_value_bounded, log_expected < floor_log)
  numeric_expected = exp(log_expected)
  numeric_expected[log_expected < floor_log] = .Machine$double.xmin
  expect_equal(result$p_value, numeric_expected, tolerance = 0)
  expect_identical(result$p_value[result$p_value_bounded], rep(.Machine$double.xmin, sum(result$p_value_bounded)))
  expect_true(all(result$p_value > 0))
  expect_true(all(diff(result$p_value) <= 0))
  expect_identical(result$inference_status, rep("testable", length(beta)))
})

test_that("constant and invalid nulls are untestable with explicit statuses", {
  result = tkoi_normal_inference(
    observed = rep(1, 7),
    null_mean = c(0, 0, 0, 0, 0, Inf, NA_real_),
    null_sd = c(0, -1, Inf, NA_real_, NaN, 1, 1)
  )
  expect_identical(result$inference_status, c("constant_null", rep("invalid_null", 6)))
  expect_true(all(is.nan(result$beta)))
  expect_true(all(is.nan(result$p_value)))
  expect_true(all(is.nan(result$log_p_value)))
  expect_true(all(is.na(result$p_value_bounded)))
})

test_that("nonfinite observed values and overflowing beta cannot become significant", {
  result = tkoi_normal_inference(
    observed = c(NA_real_, NaN, Inf, -Inf, .Machine$double.xmax),
    null_mean = c(0, 0, 0, 0, -.Machine$double.xmax),
    null_sd = rep(1, 5)
  )
  expect_identical(result$inference_status, c(rep("nonfinite_observed", 4), "nonfinite_beta"))
  expect_true(all(is.nan(result$beta)))
  expect_true(all(is.nan(result$p_value)))
  expect_true(all(is.nan(result$log_p_value)))
  expect_true(all(is.na(result$p_value_bounded)))
})

test_that("log-space BH matches p.adjust on ordinary probabilities with missing values", {
  p = withr::with_seed(812, stats::runif(100))
  p[c(2, 9, 37)] = NA_real_
  p[c(5, 15)] = c(0, 1)
  names(p) = paste0("node", seq_along(p))
  for (n in c(length(p), length(p) + 23L)) {
    expected = stats::p.adjust(p, method = "BH", n = n)
    actual = .tkoi_bh_log(log(p), n = n)
    expect_equal(actual, log(expected), tolerance = 1e-14)
    expect_identical(is.na(actual), is.na(p))
  }
  expect_equal(.tkoi_bh_log(log(p)), log(stats::p.adjust(p, "BH")), tolerance = 1e-14)
  expect_equal(.tkoi_bh_log(log(c(0.01, NA, 0.03))), log(c(0.02, NA, 0.03)), tolerance = 1e-14)
  expect_equal(.tkoi_bh_log(log(c(0.01, NA, 0.03)), n = 3), log(c(0.03, NA, 0.045)), tolerance = 1e-14)
})

test_that("log-space BH handles empty, single, and all-missing families", {
  expect_identical(.tkoi_bh_log(numeric()), numeric())
  expect_equal(.tkoi_bh_log(log(0.04)), log(0.04), tolerance = 0)
  missing = c(a = NA_real_, b = NA_real_)
  expect_identical(.tkoi_bh_log(missing), log(stats::p.adjust(missing, "BH")))
  expect_identical(.tkoi_bh_log(missing, n = 5), log(stats::p.adjust(missing, "BH", n = 5)))
})

test_that("log-space BH preserves ordering and multiplicity below double probability range", {
  log_p = c(-1000, -1001, -1002, NA_real_, -2000)
  expected = c(-1000, -1001 + log(4 / 3), -1002 + log(2), NA_real_, -2000 + log(4))
  actual = .tkoi_bh_log(log_p)
  expect_equal(actual, expected, tolerance = 1e-13)
  expect_true(all(is.finite(actual[!is.na(log_p)])))
  expect_identical(order(actual, na.last = TRUE), order(log_p, na.last = TRUE))
  expect_true(all(actual[!is.na(log_p)] >= log_p[!is.na(log_p)]))
  expect_equal(.tkoi_bh_log(c(-1000, -1000, -1002)), c(-1000, -1000, -1002 + log(3)), tolerance = 1e-13)
  expect_equal(.tkoi_bh_log(log_p, n = 8), expected + log(2), tolerance = 1e-13)
})

test_that("log-space BH never lowers an unadjusted log probability through rounding", {
  # Subtracting log(rank) after adding log(n) used to lower the last value by
  # 5.55e-17 for this family, even though its BH multiplier is exactly one.
  witness = log(seq(0.001, 0.9, length.out = 3))
  adjusted = .tkoi_bh_log(witness)
  expect_true(all(adjusted >= witness))
  expect_identical(adjusted[3], witness[3])

  families = list(
    c(-1e-30, -1e-17, -.Machine$double.eps, 0),
    log(seq(0.0001, 0.99, length.out = 101)),
    c(-1000, -1000, -1001, -2000, NA_real_)
  )
  for (log_p in families) {
    for (n in c(length(log_p), length(log_p) + 3L)) {
      adjusted = .tkoi_bh_log(log_p, n = n)
      available = !is.na(log_p)
      expect_true(all(adjusted[available] >= log_p[available]))
      expect_true(all(adjusted[available] <= 0))
      expect_identical(is.na(adjusted), is.na(log_p))
    }
  }
})

test_that("bounded probabilities preserve missing values and never replace small probabilities by zero", {
  floor_log = log(.Machine$double.xmin)
  log_p = c(0, log(0.05), floor_log, floor_log - 1, -1000, -Inf, NA_real_, NaN)
  actual = .tkoi_bounded_probability(log_p)
  numeric_expected = exp(log_p[1:6])
  numeric_expected[log_p[1:6] < floor_log] = .Machine$double.xmin
  expect_equal(actual[1:6], numeric_expected, tolerance = 0)
  expect_identical(actual[4:6], rep(.Machine$double.xmin, 3))
  expect_true(all(actual[1:6] > 0 & actual[1:6] <= 1))
  expect_true(all(is.na(actual[7:8])))
})

test_that("probability formatting distinguishes bounded values from actual floor values", {
  floor_log = log(.Machine$double.xmin)
  x = .tkoi_bounded_probability(c(log(1.2345e-5), floor_log, -1000))
  log_p = c(log(1.2345e-5), floor_log, -1000)
  plain = tkoi_format_probability(x, log_p = log_p, format = "plain", digits = 3)
  expect_identical(plain[1], "1.23e-05")
  expect_false(grepl("^\\s*<=", plain[2]))
  expect_match(plain[3], "^\\s*<=\\s*2\\.23e-308$")
  flags = tkoi_format_probability(x, format = "plain", digits = 3, bounded = c(FALSE, FALSE, TRUE))
  expect_identical(flags, plain)
})

test_that("HTML probability formatting supplies formatted exponents and bounded markers", {
  log_p = c(log(1.2345e-5), -1000)
  x = .tkoi_bounded_probability(log_p)
  html = tkoi_format_probability(x, log_p = log_p, format = "html", digits = 3)
  expect_match(html[1], "1\\.23")
  expect_match(html[1], "(×|&times;)\\s*10\\s*<sup>[-−]0?5</sup>")
  expect_match(html[2], "(≤|&le;|&leq;|&lt;=)")
  expect_match(html[2], "(×|&times;)\\s*10\\s*<sup>[-−]308</sup>")
})

test_that("stored zero probabilities are recovered only when their log values provide evidence", {
  x = c(recoverable = 0, below_floor = 0, negative_infinite_log = 0, unknown = 0, missing = NA_real_)
  log_p = c(log(1e-12), -1000, -Inf, NA_real_, NA_real_)
  original_x = x
  original_log = log_p
  plain = tkoi_format_probability(x, log_p = log_p, format = "plain", digits = 3)

  expect_identical(unname(plain), c(
    "1.00e-12", "<=2.23e-308", "<=2.23e-308",
    "Unresolved zero (log unavailable)", "NA"
  ))
  expect_identical(names(plain), names(x))
  expect_identical(x, original_x)
  expect_identical(log_p, original_log)
  expect_identical(tkoi_format_probability(0), "Unresolved zero (log unavailable)")
  expect_identical(tkoi_format_probability(0, log_p = NaN), "Unresolved zero (log unavailable)")

  html = tkoi_format_probability(x, log_p = log_p, format = "html", digits = 3)
  expect_match(unname(html[1]), "^1\\.00 (×|&times;) 10<sup>[-−]12</sup>$")
  expect_match(unname(html[2]), "^(≤|&le;|&leq;|&lt;=)2\\.23 (×|&times;) 10<sup>[-−]308</sup>$")
  expect_identical(unname(html[3]), unname(html[2]))
  expect_identical(unname(html[4]), "Unresolved zero (log unavailable)")
  expect_identical(unname(html[5]), "NA")
  expect_identical(x, original_x)
})

test_that("display tables preserve stored zeros while showing recovered or unresolved probabilities", {
  source = data.frame(p_value = c(0, 0, 0), log_p_value = c(log(1e-12), -1000, NA_real_))
  displayed = tkoi_probability_table(source)
  expect_identical(displayed$p_value, source$p_value)
  expect_identical(displayed$log_p_value, source$log_p_value)
  expect_identical(displayed$p_value_display, c(
    "1.00e-12", "<=2.23e-308", "Unresolved zero (log unavailable)"
  ))
  expect_identical(names(source), c("p_value", "log_p_value"))
})

test_that("probability formatting neither changes nor depends on global print options", {
  before = options()[c("digits", "scipen", "OutDec")]
  x = c(0.123456, 1.2345e-8, .Machine$double.xmin)
  a = withr::with_options(list(digits = 3, scipen = 999), tkoi_format_probability(x, digits = 3))
  b = withr::with_options(list(digits = 15, scipen = -9), tkoi_format_probability(x, digits = 3))
  expect_identical(a, b)
  expect_identical(options()[c("digits", "scipen", "OutDec")], before)
})
