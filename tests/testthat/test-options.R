# Tests for the optional behaviours added in this branch

outcome <- function(d, m, ...) {
  assess_replication_outcome(d, metric = m, ...)$results[[paste0(m, "_outcome")]]
}

test_that("orig_p_value_inexact = 'bound' keeps 'less-than' p-values", {
  args <- list(r_o = 0.046, n_o = 236, r_r = 0.2, n_r = 1000,
               orig_p_value = 0.001, orig_p_value_type = "less-than")
  drop  <- do.call(prepare_effect_sizes, args)
  bound <- do.call(prepare_effect_sizes, c(args, orig_p_value_inexact = "bound"))
  expect_true(is.na(drop$orig_p_value))
  expect_equal(bound$orig_p_value, 0.001)
  expect_identical(outcome(drop, "significance_r"), "OS not significant")
  expect_identical(outcome(bound, "significance_r"), "Replication success")
  gt <- do.call(prepare_effect_sizes, modifyList(args, list(orig_p_value_type = "greater-than",
                                                            orig_p_value_inexact = "bound")))
  expect_true(is.na(gt$orig_p_value))
})

test_that("repli_p_value_tails doubles one-tailed replication p-values", {
  d <- prepare_effect_sizes(r_o = 0.3, n_o = 100, r_r = 0.2, n_r = 100,
                            repli_p_value = c(0.02), repli_p_value_tails = "one-tailed")
  expect_equal(d$repli_p_value, 0.04)
})

test_that("orient_by_same_direction makes all metrics use the same direction", {
  d <- prepare_effect_sizes(r_o = c(0.3, 0.3), n_o = c(200, 200), r_r = c(-0.3, 0.3),
                            n_r = c(200, 200), same_direction = c(TRUE, FALSE),
                            orient_by_same_direction = TRUE)
  expect_equal(sign(d$es_r), c(1, -1))
  ms <- c("significance_r", "BF_JZS", "consistency_pi", "equivalence", "small_telescopes")
  o <- sapply(ms, function(m) outcome(d, m))
  expect_true(all(o[1, ] == "Replication success"))
  expect_true(all(o[2, ] == "Replication failure"))
  expect_error(prepare_effect_sizes(r_o = 0.3, n_o = 200, r_r = 0.3, n_r = 200,
                                    orient_by_same_direction = TRUE))
})

test_that("use_pattern controls the pattern flag", {
  d <- prepare_effect_sizes(r_o = c(0.3, 0.3, 0.3), n_o = rep(200, 3), r_r = c(0.3, 0.3, 0.3),
                            n_r = rep(200, 3), pattern_met = c(TRUE, FALSE, NA))
  expect_identical(outcome(d, "significance_r"),
                   c("Replication success", "Replication failure", "Replication success"))
  expect_identical(outcome(d, "consistency_pi"), rep("Replication success", 3))
  expect_identical(outcome(d, "significance_r", use_pattern = "none"), rep("Replication success", 3))
  expect_identical(outcome(d, "consistency_pi", use_pattern = "all"),
                   c("Replication success", "Replication failure", "Replication success"))
  expect_false(is.na(outcome(d, "BFr")[3]))   # missing flag no longer gives NA
})

test_that("inverse-variance meta-analysis differs from Stouffer when v_o != v_r", {
  so <- 1 / sqrt(47); sr <- 1 / sqrt(797)
  d <- prepare_effect_sizes(z_o = 4.2 * so, z_r = 0.5 * sr, se_o = so, se_r = sr,
                            n_o = 50, n_r = 800, effectsize = "ZCOR")
  expect_identical(outcome(d, "significance_agg"), "Replication success")
  expect_identical(outcome(d, "significance_agg",
                           metric_args = list(significance_agg = list(method = "iv"))),
                   "Replication failure")
})

test_that("d_to_r uses a = 4 for equal groups", {
  expect_equal(d_to_r(0.5, n1 = 20, n2 = 20)[["r"]], 0.5 / sqrt(0.25 + 4))
  expect_equal(d_to_r(0.5, n1 = 20, n2 = 20)[["r"]], d_to_r(0.5)[["r"]])
})

test_that("the default Fisher-z SE is 1/sqrt(n - 3)", {
  d <- prepare_effect_sizes(r_o = 0.5, n_o = 20, r_r = 0.5, n_r = 50)
  expect_equal(c(d$se_o, d$se_r), 1 / sqrt(c(17, 47)))
})

test_that("scalar p-value metadata applies to every pair", {
  d <- prepare_effect_sizes(r_o = c(0.3, 0.3), n_o = c(100, 100), r_r = c(0.2, 0.2),
                            n_r = c(100, 100), repli_p_value = c(0.01, 0.2),
                            repli_p_value_tails = "one-tailed")
  expect_equal(d$repli_p_value, c(0.02, 0.4))
  b <- prepare_effect_sizes(r_o = c(0.3, 0.3), n_o = c(100, 100), r_r = c(0.2, 0.2),
                            n_r = c(100, 100), orig_p_value = c(0.001, 0.01),
                            orig_p_value_type = "less-than", orig_p_value_inexact = "bound")
  expect_equal(b$orig_p_value, c(0.001, 0.01))
})

test_that("orient_by_same_direction leaves missing and zero estimates alone", {
  d <- prepare_effect_sizes(r_o = c(0.3, 0.3, 0), n_o = rep(100, 3), r_r = c(-0.3, NA, 0.2),
                            n_r = rep(100, 3), same_direction = c(TRUE, TRUE, FALSE),
                            orient_by_same_direction = TRUE)
  expect_equal(d$orientation_mirrored, c(TRUE, FALSE, FALSE))
  expect_true(is.na(d$es_r[2]))
  expect_equal(d$es_r[3], atanh(0.2))
})

test_that("inverse-variance meta-analysis uses same_direction", {
  d <- prepare_effect_sizes(z_o = 0.5, z_r = 0.5, se_o = 0.1, se_r = 0.1, n_o = 103, n_r = 103,
                            effectsize = "ZCOR", same_direction = FALSE)
  expect_identical(outcome(d, "significance_agg",
                           metric_args = list(significance_agg = list(method = "iv"))),
                   "Replication failure")
})

test_that("inverse-variance meta-analysis has overall T1E at most alpha^2 for a pre-specified direction", {
  set.seed(1)
  N <- 4e5; se_o <- 0.2; se_r <- 0.1
  zo <- rnorm(N); zr <- rnorm(N)
  d <- data.frame(study_id = seq_len(N), effectsize = "ZCOR", scale = "fisherz",
                  es_o = zo * se_o, se_o = se_o, es_r = zr * se_r, se_r = se_r,
                  n_o = 1 / se_o^2, n_r = 1 / se_r^2, se_src_o = "supplied", se_src_r = "supplied")
  s <- outcome(d, "significance_agg", metric_args = list(significance_agg = list(method = "iv"))) ==
    "Replication success"
  rate <- mean(s & zo > 0)            # pre-specified (positive) direction
  expect_lt(rate, 0.025^2 + 4 * sqrt(0.025^2 / N))   # at most alpha^2
  expect_gt(rate, 0.8 * 0.025^2)                     # and close to it
})

test_that("new arguments do not change positional calls", {
  f <- names(formals(prepare_effect_sizes))
  expect_identical(f[27], "stat_type_o")
  expect_identical(tail(f, 3), c("orig_p_value_inexact", "repli_p_value_tails",
                                 "orient_by_same_direction"))
})

test_that("orient_by_same_direction copes with an entirely missing estimate column", {
  for (sc in c("fisherz", "cor")) {
    d <- prepare_effect_sizes(r_o = 0.3, r_r = NA_real_, n_o = 100, n_r = 200, scale = sc,
                              same_direction = TRUE, orient_by_same_direction = TRUE)
    expect_false(d$orientation_mirrored)
  }
})
