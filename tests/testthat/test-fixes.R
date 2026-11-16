# Regression tests for the fixes in this branch and for documented values.

outcome <- function(d, m, ...) {
  assess_replication_outcome(d, metric = m, ...)$results[[paste0(m, "_outcome")]]
}
value <- function(d, m, ...) {
  assess_replication_outcome(d, metric = m, ...)$results[[paste0(m, "_value")]]
}

test_that("small telescopes is symmetric in the sign of the original", {
  pos <- prepare_effect_sizes(r_o = 0.4, n_o = 60, r_r = 0.4, n_r = 200)
  neg <- prepare_effect_sizes(r_o = -0.4, n_o = 60, r_r = -0.4, n_r = 200)
  expect_identical(outcome(neg, "small_telescopes"), outcome(pos, "small_telescopes"))
  expect_equal(value(neg, "small_telescopes"), value(pos, "small_telescopes"))
  # a precise replication in the opposite direction is not a success
  opp <- prepare_effect_sizes(r_o = -0.4, n_o = 60, r_r = 0.4, n_r = 200)
  expect_identical(outcome(opp, "small_telescopes"), "Replication failure")
})

test_that("small telescopes on the reported-limit path is symmetric too", {
  pos <- prepare_effect_sizes(r_o = 0.4, n_o = 60, r_r = 0.05, n_r = 200,
                              repli_conv_r_lb = -0.09, repli_conv_r_ub = 0.19)
  neg <- prepare_effect_sizes(r_o = -0.4, n_o = 60, r_r = -0.05, n_r = 200,
                              repli_conv_r_lb = -0.19, repli_conv_r_ub = 0.09)
  expect_identical(outcome(neg, "small_telescopes"), outcome(pos, "small_telescopes"))
})

test_that("rows without reported limits fall back to computed intervals", {
  d <- prepare_effect_sizes(r_o = c(0.4, 0.4), n_o = c(60, 60),
                            r_r = c(0.3, 0.3), n_r = c(200, 200),
                            orig_conv_r_lb = c(0.16, NA), orig_conv_r_ub = c(0.59, NA),
                            repli_conv_r_lb = c(0.17, NA), repli_conv_r_ub = c(0.42, NA))
  ref <- prepare_effect_sizes(r_o = 0.4, n_o = 60, r_r = 0.3, n_r = 200)
  for (m in c("consistency_cio", "consistency_cir", "small_telescopes")) {
    expect_false(anyNA(outcome(d, m)))
    expect_identical(outcome(d, m)[2], outcome(ref, m))
  }
})

test_that("small telescopes on scale = 'cor' works without reported limits", {
  d <- prepare_effect_sizes(r_o = c(0.4, 0.2), n_o = c(60, 80),
                            r_r = c(0.3, 0.1), n_r = c(200, 300), scale = "cor")
  expect_false(anyNA(outcome(d, "small_telescopes")))
})

test_that("SMD input on scale = 'cor' is returned as correlations", {
  args <- list(d_o = 1.6, d_r = 1.1, se_o = 0.05, se_r = 0.05,
               n1_o = 2000, n2_o = 2000, n1_r = 2000, n2_r = 2000, effectsize = "SMD")
  a <- do.call(prepare_effect_sizes, c(args, scale = "cor"))
  b <- do.call(prepare_effect_sizes, c(args, scale = "fisherz"))
  expect_equal(atanh(a$es_o), b$es_o)
  expect_identical(outcome(a, "equivalence"), outcome(b, "equivalence"))
})

test_that("snapshot hybrid uses the prior arguments", {
  d <- prepare_effect_sizes(r_o = 0.35, n_o = 80, r_r = 0.15, n_r = 200)
  v0 <- value(d, "snaphybrid")
  v1 <- value(d, "snaphybrid",
              metric_args = list(snaphybrid = list(pr_0 = 0.97, pr_s = 0.01, pr_m = 0.01)))
  expect_false(isTRUE(all.equal(v0, v1)))
  expect_equal(value(d, "snaphybrid",
                     metric_args = list(snaphybrid = list(pr_0 = 0.25, pr_s = 0.25, pr_m = 0.25))), v0)
  expect_error(assess_replication_outcome(
    d, metric = "snaphybrid",
    metric_args = list(snaphybrid = list(pr_0 = 0.6, pr_s = 0.3, pr_m = 0.3))))
})

test_that("sceptical p-value uses the SE ratio for supplied SEs", {
  d1 <- prepare_effect_sizes(z_o = 0.235, z_r = 0.21, se_o = 0.1, se_r = 0.1,
                             n_o = 100, n_r = 100, effectsize = "ZCOR")
  d2 <- prepare_effect_sizes(z_o = 0.235, z_r = 0.21, se_o = 0.1, se_r = 0.1,
                             n_o = 100, n_r = 400, effectsize = "ZCOR")
  expect_equal(value(d1, "pScept"), value(d2, "pScept"))
})

test_that("calibrated thresholds match the documented values", {
  expect_equal(calibrate_threshold("significance_agg", 0.025), sqrt(2) * qnorm(1 - 0.025^2))
  expect_equal(calibrate_threshold("sum_p", 0.025), sqrt(2) * 0.025)
  expect_equal(round(calibrate_threshold("BF_JZS", 0.025), 2), 1.60)
  expect_equal(round(calibrate_threshold("BF_meta", 0.025), 1), 23.7)
  expect_equal(round(calibrate_threshold("BFs", c(0.01, 0.025, 0.05)), 2), c(3.71, 2.36, 1.65))
  expect_equal(round(calibrate_threshold("BFs", 0.025, c = c(0.25, 1, 4)), 2), c(2.95, 2.36, 1.42))
})

test_that("BFe detailed labels agree with the outcome", {
  d <- prepare_effect_sizes(z_o = 0.3, z_r = 0.25, se_o = 0.1, se_r = 0.1,
                            n_o = 103, n_r = 103, effectsize = "ZCOR")
  r <- assess_replication_outcome(d, metric = "BFe", metric_args = list(BFe = list(gamma = 3)))$results
  expect_identical(r$BFe_outcome, "Inconclusive")
  expect_identical(r$BFe_outcome_detailed, "Inconclusive evidence")
})

test_that("a missing original estimate gives NA for small telescopes", {
  d <- prepare_effect_sizes(r_o = c(0.4, NA), n_o = c(60, 60), r_r = c(0.3, 0.3),
                            n_r = c(200, 200), repli_conv_r_lb = c(0.17, 0.17),
                            repli_conv_r_ub = c(0.42, 0.42))
  expect_identical(outcome(d, "small_telescopes"), c("Replication success", NA))
})

test_that("rows with missing reported limits match the computed path exactly", {
  a <- prepare_effect_sizes(r_o = 0.4, n_o = 60, r_r = 0.0615, n_r = 200)
  b <- prepare_effect_sizes(r_o = 0.4, n_o = 60, r_r = 0.0615, n_r = 200,
                            repli_conv_r_lb = NA_real_, repli_conv_r_ub = NA_real_)
  expect_identical(outcome(b, "small_telescopes"), outcome(a, "small_telescopes"))
  expect_equal(value(b, "small_telescopes"), value(a, "small_telescopes"))
})

test_that("a single reported CI side is used", {
  for (sc in c("fisherz", "cor")) {
    one  <- prepare_effect_sizes(r_o = 0.4, n_o = 60, r_r = 0, n_r = 200,
                                 repli_conv_r_ub = 0.8, scale = sc)
    both <- prepare_effect_sizes(r_o = 0.4, n_o = 60, r_r = 0, n_r = 200,
                                 repli_conv_r_ub = 0.8, repli_conv_r_lb = NA_real_, scale = sc)
    for (m in c("small_telescopes", "consistency_cir"))
      expect_identical(outcome(one, m), outcome(both, m))
    expect_identical(outcome(one, "consistency_cir"), "Replication success")
  }
})

test_that("snapshot hybrid accepts zero priors", {
  d <- prepare_effect_sizes(r_o = 0.35, n_o = 80, r_r = 0.15, n_r = 200)
  args <- function(p) list(snaphybrid = list(pr_0 = p[1], pr_s = p[2], pr_m = p[3]))
  expect_false(is.na(value(d, "snaphybrid", metric_args = args(c(0, 0.25, 0.25)))))
  expect_identical(outcome(d, "snaphybrid", metric_args = args(c(1, 0, 0))), "Replication failure")
})

test_that("length-one metadata is recycled and wrong lengths are rejected", {
  d <- prepare_effect_sizes(r_o = c(0.3, 0.4), n_o = c(100, 100), r_r = c(0.2, 0.1),
                            n_r = c(100, 100), orig_p_value = c(0.001, 0.01),
                            orig_p_value_type = "exact", orig_p_value_tails = "one-tailed")
  expect_equal(d$orig_p_value, c(0.002, 0.02))
  expect_error(prepare_effect_sizes(r_o = c(0.3, 0.4, 0.5), n_o = rep(100, 3),
                                    r_r = c(0.2, 0.1, 0.1), n_r = rep(100, 3),
                                    stat_type_o = c("t", "t")))
})

test_that("snapshot priors that sum to one up to rounding are accepted", {
  d <- prepare_effect_sizes(r_o = 0.35, n_o = 80, r_r = 0.15, n_r = 200)
  expect_false(is.na(value(d, "snaphybrid",
    metric_args = list(snaphybrid = list(pr_0 = 0.33, pr_s = 0.56, pr_m = 0.11)))))
})

test_that("reported CI vectors of the wrong length are rejected", {
  expect_error(prepare_effect_sizes(r_o = rep(0.3, 4), n_o = rep(100, 4), r_r = rep(0.2, 4),
                                    n_r = rep(200, 4), repli_conv_r_ub = c(0.4, 0.6)))
})
