# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------

# %||% is only in base R from 4.4.0; define a fallback for older versions
if (!exists("%||%")) `%||%` <- function(x, y) if (is.null(x)) y else x


#' Convert Cohen's d to Point-Biserial correlation
#'
#' Converts a standardised mean difference (Cohen's d) to a point-biserial
#' correlation coefficient using the formula from Cooper et al. (2009). When group
#' sizes \code{n1} and \code{n2} are provided, the exact correction factor is
#' used. Otherwise a balanced-groups approximation (a = 4) is assumed.
#'
#' @param d Cohen's d effect size.
#' @param se_d Standard error of Cohen's d. If omitted, only the
#'   converted effect size is returned.
#' @param n1 Sample size of group 1 (treatment group). Optional; used together with
#'   \code{n2} to compute the exact correction factor \code{a}.
#' @param n2 Sample size of group 2 (placebo group). Optional; used together with
#'   \code{n1} to compute the exact correction factor \code{a}.
#'
#' @returns A named numeric vector. If \code{se_d} is provided, returns
#'   \code{c("r", "Std. Err.")} containing the converted correlation and its
#'   delta-method standard error. If \code{se_d} is omitted, returns only
#'   \code{"r"}.
#'
#' @references
#' Cooper, H., Hedges, L. V., & Valentine, J. C. (Eds.). (2009).
#' \emph{The Handbook of Research Synthesis and Meta-Analysis} (2nd ed.).
#' Russell Sage Foundation.
#'
#' @export
#'
#' @examples
#' # Without standard error, balanced groups assumed
#' d_to_r(d = 0.5)
#'
#' # With standard error and known group sizes
#' d_to_r(d = 0.5, se_d = 0.1, n1 = 50, n2 = 50)
d_to_r <- function(d, se_d, n1=NULL, n2=NULL){
  # Handle cases where r (correlation) is not specified in the right format
  stopifnot(!missing(d), is.numeric(d))

  if (!missing(n1) & !missing(n2)){
    a <- (n1 + n2 - 2)^2 / ((n1)*(n2))
  } else{
    a <- 4
  }
  # Convert Cohen's d into correlation r
  r_d <- d / sqrt(d^2 + a)

  if (missing(se_d)){
    out <- r_d
    names(out) <- c("r")
  } else{
    r_d_se <- sqrt( (a^2 * se_d^2) / (d^2 + a)^3)

    out <- c(r_d, r_d_se)
    names(out) <- c("r", "Std. Err.")

  }
  return(out)
}

#' Convert Correlations to Fisher Z Transformation
#'
#' Applies the Fisher Z transformation (inverse-hyperbolic tangent) to a
#' point-biserial correlation, and propagates its standard error
#' using the delta method. This function is intended for correlations that
#' were previously converted from Cohen's d via \code{\link{d_to_r}}, not
#' for directly observed correlations.
#'
#' @param r Correlation coefficient to be transformed. Must be in
#'   the interval \code{[-1, 1]}.
#' @param se_r Standard error of the correlation \code{r}. Optional:
#'   if provided, the delta-method standard error of the Fisher Z transformed
#'   correlation is also returned.
#'
#' @returns A named numeric vector. If \code{se_r} is provided, returns
#'   \code{c("Fisher Z correlation", "Std. Err.")} containing the transformed
#'   correlation and its delta-method standard error. If \code{se_r} is
#'   omitted, returns only \code{"Fisher Z correlation"}.
#'
#' @details
#' The delta method approximation uses the first derivative of \code{atanh(r)},
#' which is \code{1 / (1 - r^2)}, to propagate the standard error:
#' \deqn{SE(Z_r) = SE(r) \cdot \frac{1}{1 - r^2}}
#'
#' @seealso \code{\link{d_to_r}} for converting Cohen's d to a correlation
#'   prior to applying this function.
#'
#' @export
#'
#' @examples
#' # Without standard error
#' r_to_zr(r = 0.3)
#'
#' # With standard error
#' r_to_zr(r = 0.3, se_r = 0.05)
r_to_zr <- function(r, se_r){
  # Handle cases where the correlation r is not specified in the right format
  stopifnot(!missing(r), is.numeric(r), (r >= -1 & r <=1))

  # Apply Fisher Z transf. to correlation
  zr <- atanh(r)

  if (missing(se_r)){
    out <- zr
    names(out) <- c("Fisher Z correlation")
  } else {
    # Apply delta method for fisher z corr. (1 / (1-x^2) is 1st deriv. of artanh(x))
    se_zr <- sqrt( se_r^2 * (1/(1-r^2))^2 )

    out <- c(zr, se_zr)
    names(out) <- c("Fisher Z correlation", "Std. Err.")
  }
  return(out)
}

#' Convert Log Odds Ratio to Cohen's d
#'
#' Converts a log odds ratio to a standardised mean difference (Cohen's d)
#' using the formula from Cooper et al. (2009). Optionally propagates the standard
#' error.
#'
#' @param logor Log odds ratio to be converted.
#' @param se_logor Standard error of the log odds ratio. Optional:
#'   if provided, the standard error of the converted Cohen's d is also
#'   returned.
#'
#' @returns A named numeric vector. If \code{se_logor} is provided, returns
#'   \code{c("Cohen's d", "Std. Err.")} containing the converted effect size
#'   and its standard error. If \code{se_logor} is omitted, returns only
#'   \code{"Cohen's d"}.
#'
#' @details
#' The conversion formula is:
#' \deqn{d = \log(OR) \cdot \frac{\sqrt{3}}{\pi}}
#' and the standard error propagates as:
#' \deqn{SE(d) = SE(\log(OR)) \cdot \frac{\sqrt{3}}{\pi}}
#'
#' @references
#' Cooper, H., Hedges, L. V., & Valentine, J. C. (Eds.). (2009).
#' \emph{The Handbook of Research Synthesis and Meta-Analysis} (2nd ed.).
#' Russell Sage Foundation.
#'
#' @seealso \code{\link{d_to_r}} for subsequently converting Cohen's d to a
#'   correlation, and \code{\link{r_to_zr}} for applying the Fisher Z
#'   transformation.
#'
#' @export
#'
#' @examples
#' # Without standard error
#' logor_to_d(logor = 0.5)
#'
#' # With standard error
#' logor_to_d(logor = 0.5, se_logor = 0.1)
logor_to_d <- function(logor, se_logor){
  d_or <- logor * sqrt(3) / pi

  if (missing(se_logor)){
    out <- d_or
    names(out) <- c("Cohen's d")
  } else{
    d_or_se <- sqrt(3) * se_logor / pi

    out <- c(d_or, d_or_se)
    names(out) <- c("Cohen's d", "Std. Err.")
  }

  return(out)
}

#' Convert Fisher Z correlation to two-sided p-value
#'
#' Computes a two-sided p-value from a Fisher Z transformed correlation and
#' its standard error. The function is vectorised and handles \code{NA} values
#' and invalid standard errors gracefully by returning \code{NA} for those
#' elements.
#'
#' @param zr Fisher Z transformed correlations.
#' @param se Standard errors of \code{zr}. Must be the same
#'   length as \code{zr}. Elements with \code{se <= 0} are treated as invalid
#'   and return \code{NA}.
#'
#' @returns A numeric vector of two-sided p-values of the same length as
#'   \code{zr}.
#'
#' @seealso \code{\link{r_to_zr}} for obtaining Fisher Z transformed
#'   correlations and their standard errors.
#'
#' @export
#'
#' @examples
#' # Single value
#' zr_to_p(zr = 0.3, se = 0.1)
#'
#' # Vectorised
#' zr_to_p(zr = c(0.3, 0.5, NA), se = c(0.1, 0.2, 0.1))
zr_to_p <- function(zr, se){
  if (!is.numeric(zr) | !is.numeric(se)) {
    stop("Both 'zr' and 'se' must be numeric")
  }
  if (length(zr) != length(se)) {
    stop("Length of 'zr' and 'se' must be the same")
  }

  p_value <- rep(NA, length(zr))

  # Handle valid cases (non-NA, se <= 0)
  valid_indices <- !(is.na(zr) | is.na(se) | se <= 0)

  if (any(valid_indices)) {
    z_value <- zr[valid_indices]
    p_value[valid_indices] <- 2 * stats::pnorm(abs(z_value), sd = se[valid_indices], lower.tail=FALSE)
  }

  return(p_value)
}

# Performs fixed effect meta-analysis on one original and one replication study
# `level`: confidence level (in %) of the pooled CI, passed to metafor::rma()
compare_effectsizes <- function(es_o, se_o, es_r, se_r, level = 95) {

  process_pair <- function(es_o, se_o, es_r, se_r) {
    if (any(is.na(c(es_o, se_o, es_r, se_r))))
      return(data.frame(
        est  = NA_real_,
        ci_lb = NA_real_,
        ci_ub = NA_real_,
        QEp  = NA_real_,
        pval = NA_real_
      ))

    # # Calculate effect sizes for each study
    # r1_result <- metafor::escalc(ni = n1, ri = r1, measure = "COR")
    # r2_result <- metafor::escalc(ni = n2, ri = r2, measure = "COR")

    # Meta-analysis on the two effect sizes
    res <- metafor::rma(yi = c(es_o, es_r),
                        sei = c(se_o, se_r),
                        method = "FE",
                        level = level)

    return(data.frame(est = res$beta[1], ci_lb = res$ci.lb, ci_ub = res$ci.ub, QEp = res$QEp, pval = res$pval))
  }

  results <- mapply(process_pair, es_o, se_o, es_r, se_r, SIMPLIFY = FALSE)

  do.call(rbind, results)

}

# Computes (95%)-CI from Fisher Z Effect Size
compute_ci <- function(es, se, alpha = 0.025){
  # Ensure es and se are the same length
  if (length(es) != length(se)) {
    stop("The length of 'es' and 'se' must be the same")
  }
  # Handle NA/non-finite values directly to avoid warnings
  valid_indices <- !is.na(es) & !is.na(se) & is.finite(es) & is.finite(se) & se > 0
  # Initialize vectors to NA
  z_lower <- z_upper <- rep(NA_real_, length(es))

  if (any(valid_indices)) {
    # Confidence interval for z (es); subset BOTH sides so lengths match
    crit <- stats::qnorm(alpha, lower.tail = FALSE)
    z_lower[valid_indices] <- es[valid_indices] - crit * se[valid_indices]
    z_upper[valid_indices] <- es[valid_indices] + crit * se[valid_indices]
  }

  # Return the confidence intervals in a data frame
  data.frame(es = es, se = se, lower = z_lower, upper = z_upper)
}

# Rescale reported 95% confidence limits to a (1 - 2 alpha) interval ----
# Reported limits cannot be recomputed, so each side is rescaled on the
# Fisher-z scale around the point estimate: distance to the limit *
# qnorm(1 - alpha) / qnorm(1 - ref_alpha). Asymmetric intervals stay
# asymmetric and a single missing limit only affects its own side.
# At alpha = ref_alpha the limits are returned unchanged.
# `on_z = TRUE`: limits and estimate are already on the Fisher-z scale.
rescale_ci <- function(lb, ub, est, alpha, on_z = FALSE, ref_alpha = 0.025) {
  if (isTRUE(all.equal(alpha, ref_alpha))) return(list(lower = lb, upper = ub))
  to_z <- function(b) atanh(pmax(pmin(b, 1 - 1e-12), -1 + 1e-12))
  zl <- if (on_z) lb  else to_z(lb)
  zu <- if (on_z) ub  else to_z(ub)
  ze <- if (on_z) est else to_z(est)
  k  <- stats::qnorm(1 - alpha) / stats::qnorm(1 - ref_alpha)
  lo <- ze - (ze - zl) * k
  hi <- ze + (zu - ze) * k
  if (on_z) list(lower = lo, upper = hi) else list(lower = tanh(lo), upper = tanh(hi))
}

# Back-transform a Fisher-z effect size + SE to the correlation scale ----
# (delta method: d/dz tanh(z) = 1 - tanh(z)^2)
fisherz_to_r <- function(z, se_z){
  r <- tanh(z)
  list(r = r, se_r = se_z * (1 - r^2))
}

# Difference test ----
diff_test <- function(es_o, es_r, se_o, se_r, alpha = 0.05){
  se    <- sqrt(se_o^2 + se_r^2)
  delta <- abs(es_r - es_o)
  p     <- stats::pnorm(delta / se, lower.tail = FALSE)   # one-tailed
  list(sig = (p <= alpha / 2), pval = 2 * p, delta = delta)
}

# Equivalence test ----
equiv_test <- function(es_o, es_r, se_o, se_r, h, alpha = 0.05){
  se    <- sqrt(se_o^2 + se_r^2)
  delta <- es_o - es_r
  z_eq  <- (delta + c(-h, h)) / se
  p_eq  <- c(stats::pnorm(z_eq[1], lower.tail = TRUE),
             stats::pnorm(z_eq[2], lower.tail = FALSE))
  list(sig = all(p_eq <= alpha), pval = max(p_eq), delta = delta)
}

# Correspondence test ----
# classification from the joint result of diff_test and equiv_test
corresp_test <- function(es_o, es_r, se_o, se_r, h, alpha = 0.05){
  d <- diff_test(es_o, es_r, se_o, se_r, alpha)
  e <- equiv_test(es_o, es_r, se_o, se_r, h, alpha)
  res <- ifelse(d$sig,
                ifelse(e$sig, "Trivial Difference", "Difference"),
                ifelse(e$sig, "Equivalence", "Indeterminacy"))
  list(ct = res, p_diff = d$pval, p_equiv = e$pval)
}


# Computes PI from Fisher Z Effect Size. Approach based on "J. M. Schauer and L. V. Hedges (2021)"
compute_prediction_interval <- function(es_o, se_o, se_r, alpha = 0.025){

  # Standard error for PI
  se_pi <- sqrt(se_o^2 + se_r^2)

  # Calculate lower and upper Bounds for PI
  PI_l <- es_o - stats::qnorm(alpha, lower.tail = FALSE) * se_pi
  PI_u <- es_o + stats::qnorm(alpha, lower.tail = FALSE) * se_pi

  return(list(lower = PI_l, upper = PI_u))
}

# Finds the effect size that the original study might have detected for a specific power level for small telescopes
find_z_for_power <- function(power, n1, n2, alpha=0.025){
  # Author of the metric suggests a 90\% confidence interval (cf. Zwaan)
  es_power <- pwr::pwr.t2n.test(power = power, n1=n1, n2=n2, sig.level = alpha,
                                alternative = "greater")
  return(es_power$d)
}


# Snapshot Hybrid. Code was taken from original authors and readapted for our needs.
hybrid_snap <- function(es, se, alpha=0.025, sm=0.1, me=0.3, la=0.5) {
  # Stop if the inputs are not a vector length 2 (original, replication)
  stopifnot(length(es) == 2, length(se) == 2)

  # Calculate critical value alpha for truncation
  ycvi <- stats::qnorm(alpha, lower.tail = FALSE) * se

  ### Transform small, medium, and large effect size to Fisher-z
  fis.sm <- r_to_zr(r = sm)[[1]]
  fis.me <- r_to_zr(r = me)[[1]]
  fis.la <- r_to_zr(r = la)[[1]]

  ### Compute probabilities
  log_f.0 <- (stats::dnorm(es[1],0,se[1], log = TRUE) -
                stats::pnorm(ycvi[1],0,se[1],lower.tail=FALSE, log.p = TRUE)) +
    stats::dnorm(es[2],0,se[2], log = TRUE)
  log_f.sm <- (stats::dnorm(es[1],fis.sm,se[1], log = TRUE) -
                 stats::pnorm(ycvi[1],fis.sm,se[1], lower.tail=FALSE, log.p = TRUE)) +
    stats::dnorm(es[2],fis.sm,se[2], log = TRUE)
  log_f.me <- (stats::dnorm(es[1],fis.me,se[1], log = TRUE) -
                 stats::pnorm(ycvi[1],fis.me,se[1], lower.tail=FALSE, log.p = TRUE)) +
    stats::dnorm(es[2],fis.me,se[2], log = TRUE)
  log_f.la <- (stats::dnorm(es[1],fis.la,se[1], log = TRUE) -
                 stats::pnorm(ycvi[1],fis.la,se[1], lower.tail=FALSE, log.p = TRUE)) +
    stats::dnorm(es[2],fis.la,se[2], log = TRUE)

  logs <- c(log_f.0, log_f.sm, log_f.me, log_f.la)
  if (any(!is.finite(logs))) {
    return(data.frame(
      p.0  = NA_real_,
      p.sm = NA_real_,
      p.me = NA_real_,
      p.la = NA_real_
    ))
  }

  # Shift by max to avoid overflow/underflow, then normalise
  m     <- max(logs)
  probs <- exp(logs - m)
  probs <- probs / sum(probs)

  return(data.frame(
    p.0  = probs[1],
    p.sm = probs[2],
    p.me = probs[3],
    p.la = probs[4]
  ))
}

#   ## Posterior probabilities of hypotheses
#
#   BF.0 <- p.0 / (1 - p.0) * 1 # Posterior odds divided by prior odds
#   BF.sm <- p.sm / (1 - p.sm) * 1 # Posterior odds divided by prior odds
#   BF.me <- p.me / (1 - p.me) * 1 # Posterior odds divided by prior odds
#   BF.la <- p.la / (1 - p.la) * 1 # Posterior odds divided by prior odds
#
#   return(data.frame(p.0 = p.0, p.sm = p.sm, p.me = p.me, p.la = p.la))
#
# }



# Calibrated success thresholds for Family 1 metrics ----

#' calibrate_threshold: Success thresholds with overall T1E rate \eqn{\alpha^2}
#'
#' Computes, for a one-sided level \eqn{\alpha}, the success threshold of a
#' Family 1 metric such that the overall T1E rate (both the original
#' and the replication effect equal to zero, the intersection null) is
#' \eqn{\alpha^2}. These are the thresholds used by
#' [assess_replication_outcome()] when `calibrated = TRUE`, so success rates
#' can be recalibrated at every value of \eqn{\alpha}.
#'
#' @param metric One of `"significance_agg"`, `"sum_p"`, `"BF_JZS"`,`"BF_meta"` or `"BFs"`.
#' @param alpha Numeric vector of one-sided levels in (0, 0.5). Default 0.025.
#' @param c Variance ratio \eqn{c = \sigma_o^2 / \sigma_r^2} used for the
#'   `"BFs"` calibration (ignored for the other metrics). A vector gives one
#'   threshold per value; [assess_replication_outcome()] uses the variance
#'   ratio of each pair. Default 1 (equal variances).
#'
#' @details
#' \describe{
#'   \item{`"significance_agg"`}{Threshold for the unweighted Stouffer sum
#'     \eqn{z_o + z_r}: \eqn{\sqrt{2} z_{1-\alpha^2}}.}
#'   \item{`"sum_p"`}{Threshold for the sum of the two one-sided p-values:
#'     \eqn{\alpha \sqrt{2}}, since \eqn{P(p_o + p_r \le s) = s^2 / 2}.}
#'   \item{`"BF_JZS"`}{Smallest threshold \eqn{\gamma} for which
#'     (\eqn{P(BF_o > \gamma)} and \eqn{P(BF_r > \gamma)}) \eqn{(\leq \alpha^2)}
#'     under the null hypothesis for every standard error, so that the overall
#'     T1E rate is at most \eqn{\alpha^2}. With the Cauchy(0, 1) prior on the
#'     effect, the critical \eqn{z} at which the Bayes factor equals
#'     \eqn{\gamma} depends on the standard error, so no single \eqn{\gamma}
#'     gives controls the overall T1E to \eqn{\alpha^2} exactly for every study.
#'     The threshold \eqn{\gamma} is found by minimising the critical \eqn{z}
#'     over the standard error and root finding on \eqn{\gamma}.}
#'   \item{`"BF_meta"`}{Smallest threshold \eqn{\gamma} for which the
#'     fixed-effect meta-analytic Bayes factor, whose success also requires
#'     the pooled estimate to have the direction of the original, gives
#'     \eqn{P(BF > \gamma) \le \alpha^2} under the null hypothesis for every
#'     standard error of the pooled estimate.}
#'   \item{`"BFs"`}{Threshold \eqn{\gamma} for the inverted sceptical Bayes
#'     factor such that \eqn{P(1/BF_S > \gamma) = \alpha^2} for a pair of
#'     studies with variance ratio \eqn{c}; pairs for which the sceptical Bayes
#'     factor does not exist count as non-successes, as in
#'     [assess_replication_outcome()]. The T1E rate depends on the data only
#'     through \eqn{c}, so a threshold computed from the pair's own \eqn{c}
#'     gives exactly \eqn{\alpha^2} for every pair, whereas a single threshold
#'     does so only at one value of \eqn{c} (e.g., \eqn{c} = 1). For a given
#'     \eqn{z_o} the success region in \eqn{z_r} is available in closed form
#'     (Pawel and Held, 2022), so the T1E is a one-dimensional
#'     integral and \eqn{\gamma} is found by root finding, without simulation.}
#' }
#' Results are cached, so repeated calls with the same arguments are fast.
#'
#' @return Numeric vector of thresholds, one per value of `alpha` (for
#'   `"BFs"`, one per combination of `alpha` and `c`, recycled to the longer
#'   of the two).
#'
#' @references
#' Pawel, S., & Held, L. (2022). The sceptical Bayes factor for the
#' assessment of replication success. \emph{Journal of the Royal Statistical
#' Society: Series B}, \bold{84}, 879-911.
#'
#' @export
#'
#' @examples
#' calibrate_threshold("BF_JZS", alpha = 0.025)             # 1.60
#' calibrate_threshold("BF_meta", alpha = 0.025)            # 23.7
#' calibrate_threshold("BFs", alpha = c(0.01, 0.025, 0.05)) # 4.38 2.81 1.98 (c = 1)
#' calibrate_threshold("BFs", alpha = 0.025, c = c(0.25, 1, 4)) # 3.57 2.81 1.58
calibrate_threshold <- function(metric = c("significance_agg", "sum_p", "BF_JZS", "BF_meta", "BFs"),
                                alpha = 0.025, c = 1) {
  metric <- match.arg(metric)
  stopifnot(is.numeric(alpha), length(alpha) >= 1, all(is.finite(alpha)),
            all(alpha > 0 & alpha < 0.5),
            is.numeric(c), length(c) >= 1, all(is.finite(c)), all(c > 0))

  key <- function(...) paste(..., sep = "|")
  switch(metric,
         significance_agg = sqrt(2) * stats::qnorm(1 - alpha^2),
         sum_p            = alpha * sqrt(2),
         BF_JZS = vapply(alpha, function(a)
           .cal_cached(key("BF_JZS", sprintf("%.17g", a)), .calibrate_bfjzs(a)), numeric(1)),
         BF_meta = vapply(alpha, function(a)
           .cal_cached(key("BF_meta", sprintf("%.17g", a)), .calibrate_bfmeta(a)), numeric(1)),
         BFs    = mapply(function(a, cc)
           .cal_cached(key("BFs", sprintf("%.17g", cc), sprintf("%.17g", a)),
                       .calibrate_bfs(a, cc)), alpha, c, USE.NAMES = FALSE)
  )
}

# internal helpers for calibrate_threshold() ----

.cal_cache <- new.env(parent = emptyenv())

# `expr` is only evaluated when the key is not cached yet
.cal_cached <- function(key, expr) {
  if (!exists(key, envir = .cal_cache, inherits = FALSE))
    assign(key, expr, envir = .cal_cache)
  get(key, envir = .cal_cache, inherits = FALSE)
}

# log JZS Bayes factor BF10 (Cauchy(0, 1) prior on the effect) for an
# estimate es with standard error se. The integrand is scaled by its value at
# es (no overflow for large z) and integrated over es +/- 40 se, where the
# likelihood has its mass; integrating over the real line can miss the peak
# when se is small.
.bfjzs_log <- function(es, se) {
  lf <- function(t) stats::dnorm(es, t, se, log = TRUE) + stats::dcauchy(t, 0, 1, log = TRUE)
  m  <- lf(es)
  I  <- stats::integrate(function(t) exp(lf(t) - m),
                         lower = es - 40 * se, upper = es + 40 * se,
                         rel.tol = 1e-10)$value
  m + log(I) - stats::dnorm(es, 0, se, log = TRUE)
}

# Critical z at which the JZS Bayes factor equals gamma (> 1) for standard
# error se. BF10 < 1 at z = 0 and increases with z, so the root is unique.
.bfjzs_crit <- function(se, gamma) {
  stats::uniroot(function(z) .bfjzs_log(z * se, se) - log(gamma),
                 lower = 0, upper = 50, tol = 1e-10)$root
}

# Smallest gamma for which the T1E rate P(BF > gamma) of a JZS Bayes factor
# test is at most `target` for every standard error (worst case over se).
#   sides = 2: success for |z| > z_crit (no direction check)
#   sides = 1: success for z > z_crit only (direction check: estimate in the
#              pre-specified direction)
.calibrate_bf_worstcase <- function(target, sides) {
  worst_rate <- function(lg) {
    zmin <- stats::optimize(function(ls) .bfjzs_crit(exp(ls), exp(lg)),
                            interval = log(c(0.01, 10)), tol = 1e-8)$objective
    sides * stats::pnorm(zmin, lower.tail = FALSE)
  }
  lo <- log(1 + 1e-6)
  if (worst_rate(lo) <= target) return(1)   # very large alpha: any gamma > 1 suffices
  exp(stats::uniroot(function(lg) worst_rate(lg) - target,
                     lower = lo, upper = log(1e8), tol = 1e-10)$root)
}

# JZS Bayes factor with a direction check: per study P(BF > gamma, estimate in
# the pre-specified direction) <= alpha, so that the overall T1E rate (both
# studies in that direction) is at most alpha^2
.calibrate_bfjzs <- function(alpha) .calibrate_bf_worstcase(alpha, sides = 1)

# Meta-analytic Bayes factor: one pooled test with a direction check,
# P(BF > gamma, pooled estimate in the claimed direction) <= alpha^2
.calibrate_bfmeta <- function(alpha) .calibrate_bf_worstcase(alpha^2, sides = 1)


# log Bayes factors of BayesRep::BFo() and BayesRep::BFr() (z-scale versions)
.bfs_logBFo <- function(zo, g) {
  stats::dnorm(zo, 0, 1, log = TRUE) - stats::dnorm(zo, 0, sqrt(1 + g), log = TRUE)
}
.bfs_logBFr <- function(zo, zr, c, g) {
  stats::dnorm(zr, 0, sqrt(1 + c * g), log = TRUE) -
    stats::dnorm(zr, zo * sqrt(c), sqrt(1 + c), log = TRUE)
}

# P(1 / BFs > gamma | zo) under H0 (zr ~ N(0, 1)), vectorised over zo.
# Success at level gamma: the original is convincing enough (minimal BFo below
# 1/gamma) and the replication BF of the sceptical prior with BFo = 1/gamma
# (closed form via the Lambert W function) against the advocacy prior is below
# 1/gamma, which is a quadratic inequality in zr.
.bfs_success_given_zo <- function(zo, gamma, c = 1) {
  out <- numeric(length(zo))
  lg  <- log(gamma)
  z2  <- zo^2
  ok  <- z2 > 1 & .bfs_logBFo(zo, pmax(z2 - 1, 0)) < -lg
  if (!any(ok)) return(out)
  zz <- z2[ok]
  K  <- zz - 2 * lg - log(zz)
  g  <- zz / (-lamW::lambertWm1(-exp(-K))) - 1     # sceptical relative prior variance
  v0 <- 1 + c * g
  v1 <- 1 + c
  m  <- sqrt(zz * c)
  A  <- 1 / v1 - 1 / v0
  B  <- -2 * m / v1
  C  <- m^2 / v1 + log(v1 / v0) + 2 * lg          # success: A zr^2 + B zr + C < 0
  p    <- numeric(length(A))
  lin  <- abs(A) < 1e-12
  disc <- B^2 - 4 * A * C
  sq   <- sqrt(pmax(disc, 0))
  r1 <- (-B - sq) / (2 * A)
  r2 <- (-B + sq) / (2 * A)
  lo <- pmin(r1, r2)
  hi <- pmax(r1, r2)
  i <- !lin & A > 0 & disc > 0
  p[i] <- stats::pnorm(lo[i], lower.tail = FALSE) - stats::pnorm(hi[i], lower.tail = FALSE)
  i <- !lin & A < 0 & disc <= 0
  p[i] <- 1
  i <- !lin & A < 0 & disc > 0
  p[i] <- stats::pnorm(lo[i]) + stats::pnorm(hi[i], lower.tail = FALSE)
  i <- lin
  p[i] <- ifelse(B[i] < 0, stats::pnorm(-C[i] / B[i], lower.tail = FALSE),
                 stats::pnorm(-C[i] / B[i]))
  out[ok] <- p
  out
}

# Overall P(1 / BFs > gamma) under H0 over all pairs, both directions
.bfs_success_prob <- function(gamma, c = 1) {
  lg   <- log(gamma)
  zmin <- stats::uniroot(function(z) .bfs_logBFo(z, z^2 - 1) + lg,
                         c(1 + 1e-9, 60), tol = 1e-12)$root
  stats::integrate(function(z) stats::dnorm(z) * .bfs_success_given_zo(z, gamma, c),
                   lower = zmin, upper = max(zmin + 1, 20), rel.tol = 1e-9)$value
}

# Root finding on log(gamma): the T1E rate decreases in gamma, each
# evaluation is one integral (a few milliseconds)
.calibrate_bfs <- function(alpha, c = 1) {
  target <- alpha^2
  f  <- function(lg) log(.bfs_success_prob(exp(lg), c)) - log(target)
  lo <- log(1.001)
  hi <- log(1e5)
  if (f(lo) < 0 || f(hi) > 0) {
    warning("No sceptical Bayes factor threshold gives the target T1E rate ",
            "for alpha = ", alpha, "; returning NA.")
    return(NA_real_)
  }
  exp(stats::uniroot(f, c(lo, hi), tol = 1e-10)$root)
}
