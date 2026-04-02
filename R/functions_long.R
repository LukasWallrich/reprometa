# Draft functions for R Package

###############################################################################
#############################  ReproMeta  #####################################
###############################################################################
library("pacman")
pacman::p_load(dplyr,FReD,ReplicationSuccess,BayesRep,checkmate,metafor)
data("RProjects")
test_df <- RProjects[RProjects$`project` == "Psychology", ]


# Main idea: Meta-function taking correlation as input primarily. 
# If correlation not available, then system recognizes cohen's d (or whatever) 
# Effect sizes will successively be converted to Fisher z correlations.

# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------

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
compare_effectsizes <- function(es_o, se_o, es_r, se_r) {
  
  process_pair <- function(es_o, se_o, es_r, se_r) {
    if (any(is.na(c(es_o, se_o, es_r, se_r)))) return(data.frame(est = NA, ci_lb = NA, QEp = NA))
    
    # # Calculate effect sizes for each study
    # r1_result <- metafor::escalc(ni = n1, ri = r1, measure = "COR")
    # r2_result <- metafor::escalc(ni = n2, ri = r2, measure = "COR")
    
    # Meta-analysis on the two effect sizes
    res <- metafor::rma(yi = c(es_o, es_r),
                        sei = c(se_o, se_r),
                        method = "FE")
    
    return(data.frame(est = res$beta[1], ci_lb = res$ci.lb, QEp = res$QEp, pval = res$pval))
  }
  
  results <- mapply(process_pair, es_o, se_o, es_r, se_r, SIMPLIFY = FALSE)
  
  do.call(rbind, results)
  
}

# Computes (95%)-CI from Fisher Z Effect Size
compute_ci <- function(es, se, alpha = 0.05){
  # Ensure es and se are the same length
  if (length(es) != length(se)) {
    stop("The length of 'es' and 'se' must be the same")
  }
  # Handle NA values directly to avoid warnings
  valid_indices <- !(is.na(es) | is.na(se) | se <= 0)
  # Initialize vectors to NA
  z_lower <- z_upper <- rep(NA, length(es))
  
  if (any(valid_indices)) {
    # Confidence interval for z (es)
    z_lower[valid_indices] <- es - stats::qnorm(alpha / 2, lower.tail = FALSE) * se
    z_upper[valid_indices] <- es + stats::qnorm(alpha / 2, lower.tail = FALSE) * se
  }
  
  # Return the confidence intervals in a data frame
  data.frame(es = es, se = se, lower = z_lower, upper = z_upper)
}

# Computes PI from Fisher Z Effect Size. Approach based on "J. M. Schauer and L. V. Hedges (2021)"
compute_prediction_interval <- function(es_o, se_o, se_r, alpha = 0.05){

  # Standard error for PI
  se_pi <- sqrt(se_o^2 + se_r^2)
  
  # Calculate lower and upper Bounds for PI
  PI_l <- es_o - stats::qnorm(alpha/2, lower.tail = FALSE) * se_pi
  PI_u <- es_o + stats::qnorm(alpha/2, lower.tail = FALSE) * se_pi
  
  return(list(lower = PI_l, upper = PI_u))
}

# Finds the effect size that the original study might have detected for a specific power level for small telescopes
find_z_for_power <- function(power, n1, n2, alpha=0.05){
  # Author of the metric suggests a 90\% confidence interval (cf. Zwaan)
  es_power <- pwr::pwr.t2n.test(power = power, n1=n1, n2=n2, sig.level = alpha)
  return(es_power$d)
}


# Snapshot Hybrid. Code was taken from original authors and readapted for our needs.
hybrid_snap <- function(es, se, alpha=0.025, sm=0.1, me=0.3, la=0.5) {
  # Stop if the inputs are not a vector length 2 (original, replication)
  stopifnot(length(es) == 2, length(se) == 2)
  
  # Calculate critical value alpha for truncation
  ycvi <- stats::qnorm(alpha, lower.tail = FALSE)*se
  
  ### Transform small, medium, and large effect size to Fisher-z
  fis.sm <- r_to_zr(r = sm)[[1]]
  fis.me <- r_to_zr(r = me)[[1]]
  fis.la <- r_to_zr(r = la)[[1]]
  
  ### Compute probabilities
  f.0 <- (stats::dnorm(es[1],0,se[1])/stats::pnorm(ycvi[1],0,se[1],
                                     lower.tail=FALSE))*stats::dnorm(es[2],0,se[2])
  f.sm <- (stats::dnorm(es[1],fis.sm,se[1])/stats::pnorm(ycvi[1],fis.sm,se[1],
                                           lower.tail=FALSE)) * stats::dnorm(es[2],fis.sm,se[2])
  f.me <- (stats::dnorm(es[1],fis.me,se[1])/stats::pnorm(ycvi[1],fis.me,se[1],
                                           lower.tail=FALSE)) * stats::dnorm(es[2],fis.me,se[2])
  f.la <- (stats::dnorm(es[1],fis.la,se[1])/stats::pnorm(ycvi[1],fis.la,se[1],
                                           lower.tail=FALSE)) * stats::dnorm(es[2],fis.la,se[2])
  
  ## Posterior probabilities of hypotheses
  p.0 <- f.0/(f.0+f.sm+f.me+f.la)
  p.sm <- f.sm/(f.0+f.sm+f.me+f.la)
  p.me <- f.me/(f.0+f.sm+f.me+f.la)
  p.la <- f.la/(f.0+f.sm+f.me+f.la)
  
  BF.0 <- p.0 / (1 - p.0) * 1 # Posterior odds divided by prior odds
  BF.sm <- p.sm / (1 - p.sm) * 1 # Posterior odds divided by prior odds
  BF.me <- p.me / (1 - p.me) * 1 # Posterior odds divided by prior odds
  BF.la <- p.la / (1 - p.la) * 1 # Posterior odds divided by prior odds
  
  return(data.frame(p.0 = p.0, p.sm = p.sm, p.me = p.me, p.la = p.la))
  
}


#####################
### Meta-function ### ----
#####################

#' Assess Replication outcomes using multiple statistical metrics
#'
#' Evaluates the outcome of replication studies using one or more statistical
#' metrics. Effect sizes from original and replication studies are first
#' converted to Fisher Z transformed correlations, then passed to the selected
#' metrics. Results are returned in a long-format data frame with one row per
#' study per metric.
#' 
#' @param r_o Correlation effect size(s) of the original studies.
#'   Required when \code{effectsize = "COR"} (Correlation is set as default).
#' @param n_o Sample size(s) of the original studies. Required
#'   when \code{effectsize = "COR"}.
#' @param r_r Correlation effect size(s) of the replication studies. 
#'   Required when \code{effectsize = "COR"} (Correlation is set as default).
#' @param n_r Sample size(s) of the replication studies. Required
#'   when \code{effectsize = "COR"}.
#' @param metric One or more statistical metrics to use for assessing
#'   replication outcomes. Available metrics are:
#'   \describe{
#'     \item{\code{"significance_r"}}{Significance criterion (vote counting):
#'       replication is a success if the replication effect is statistically
#'       significant and in the same direction as the original.}
#'     \item{\code{"significance_agg"}}{Fixed-effects meta-analysis aggregating
#'       the original and replication effect sizes.}
#'     \item{\code{"consistency_cir"}}{Checks whether the original effect size
#'       falls within the replication's confidence interval.}
#'     \item{\code{"consistency_cio"}}{Checks whether the replication effect
#'       size falls within the original's confidence interval.}
#'     \item{\code{"consistency_pi"}}{Checks whether the replication effect
#'       size falls within the original's prediction interval.}
#'     \item{\code{"small_telescopes"}}{Small telescopes approach: replication
#'       is a success if the original study had at least 33\% power to detect
#'       the replication effect.}
#'     \item{\code{"BF_JZS"}}{Jeffreys-Zellner-Siow (JZS) Bayes factor
#'       comparing evidence for a non-zero effect against the null in both
#'       studies.}
#'     \item{\code{"BFe"}}{Equality of effects Bayes factor testing whether
#'       original and replication effects are equal.}
#'     \item{\code{"BFr"}}{Replication Bayes factor using the original effect
#'       as the alternative hypothesis for the replication.}
#'     \item{\code{"BFs"}}{Sceptical Bayes factor testing replication success
#'       against a sceptical prior. Returns \code{NaN} when evidence is too
#'       inconclusive to establish replication success at any level.}
#'     \item{\code{"pScept"}}{Controlled sceptical p-value; replication is a
#'       success if \eqn{p_s \leq \alpha}.}
#'     \item{\code{"snaphybrid"}}{Snapshot hybrid method; classifies the
#'       replication into null, small, medium, or large effect based on
#'       posterior probabilities.}
#'   }
#' @param effectsize The type of effect size provided.
#'   One of \code{"COR"} (Pearson correlation, default), \code{"SMD"}
#'   (standardised mean difference / Cohen's d), \code{"LOGOR"} (log odds
#'   ratio), or \code{"ZCOR"} (Fisher Z transformed correlation).
#' @param d_o Cohen's d effect size(s) of the original study(s).
#'   Required when \code{effectsize = "SMD"}.
#' @param d_r Cohen's d effect size(s) of the replication study(s).
#'   Required when \code{effectsize = "SMD"}.
#' @param se_o Standard error(s) of the original effect size(s).
#'   Required when \code{effectsize} is \code{"SMD"}, \code{"LOGOR"}, or
#'   \code{"ZCOR"}.
#' @param se_r Standard error(s) of the replication effect size(s).
#'   Required when \code{effectsize} is \code{"SMD"}, \code{"LOGOR"}, or
#'   \code{"ZCOR"}.
#' @param logor_o Log odds ratio(s) of the original study(s).
#'   Required when \code{effectsize = "LOGOR"}.
#' @param logor_r Log odds ratio(s) of the replication study(s).
#'   Required when \code{effectsize = "LOGOR"}.
#' @param n1_o Sample size(s) of group 1 (treatment) in the original study(s).
#'   Required when \code{effectsize} is \code{"SMD"} or \code{"LOGOR"}.
#' @param n2_o Sample size(s) of group 2 (placebo) in the original study(s).
#'   Required when \code{effectsize} is \code{"SMD"} or \code{"LOGOR"}.
#' @param n1_r Sample size(s) of group 1 (treatment) in the replication study(s).
#'   Required when \code{effectsize} is \code{"SMD"} or \code{"LOGOR"}.
#' @param n2_r Sample size(s) of group 2 (placebo) in the replication study(s).
#'   Required when \code{effectsize} is \code{"SMD"} or \code{"LOGOR"}.
#' @param z_o Fisher Z transformed correlation(s) of the original
#'   study(s). Required when \code{effectsize = "ZCOR"}.
#' @param z_r Fisher Z transformed correlation(s) of the
#'   replication study(s). Required when \code{effectsize = "ZCOR"}.
#'
#' @returns A data frame in long format with one row per study per metric,
#'   ordered by \code{study_id} and \code{metric}, containing the following
#'   columns:
#'   \describe{
#'     \item{\code{study_id}}{Integer. Row index of the study in the input.}
#'     \item{\code{metric}}{Character. Name of the metric used.}
#'     \item{\code{quantitative_value}}{Numeric. The quantitative output of the
#'       metric (e.g. p-value, Bayes factor, posterior probability).}
#'     \item{\code{es_o}}{Numeric. Fisher Z transformed effect size of the
#'       original study.}
#'     \item{\code{se_o}}{Numeric. Standard error of the original Fisher Z
#'       effect size.}
#'     \item{\code{es_r}}{Numeric. Fisher Z transformed effect size of the
#'       replication study.}
#'     \item{\code{se_r}}{Numeric. Standard error of the replication Fisher Z
#'       effect size.}
#'     \item{\code{outcome}}{Character. Replication outcome classification
#'       (e.g. \code{"Replication success"}, \code{"Replication failure"},
#'       \code{"Inconclusive"}).}
#'     \item{\code{outcome_detailed}}{Character. A more detailed description
#'       of the outcome for the given metric.}
#'   }
#' 
#' @seealso
#'   \code{\link{d_to_r}}, \code{\link{r_to_zr}}, \code{\link{logor_to_d}},
#'   \code{\link{zr_to_p}}
#'   
#' @export
#'
#' @examples
#' # Single metric, correlation input
#' assess_replication_outcome(
#'   d_o = 0.35, se_o = 0.082, n1o = 42, n2o = 40,
#'   d_r = 0.26, se_r = 0.110, n1r = 30, n2r = 30,
#'   effectsize = "SMD",
#'   metric = "significance_r"
#' )
#'
#' # Multiple metrics
#' assess_replication_outcome(
#'   r_o = c(0.35, 0.50), n_o = c(100, 80),
#'   r_r = c(0.21, 0.10), n_r = c(120, 90),
#'   effectsize = "COR",
#'   metric = c("significance_r", "BFr", "pScept")
#' )
assess_replication_outcome <- function(r_o, n_o, r_r, n_r, metric, effectsize = "COR", 
                                       d_o=NULL, d_r=NULL, se_o=NULL, se_r=NULL,
                                       logor_o=NULL, logor_r=NULL, n1_o=NULL, n2_o=NULL,
                                       n1_r=NULL, n2_r=NULL, z_o=NULL, z_r=NULL) {
  
  # Check that relevant arguments are characters and not missing
  stopifnot(is.character(metric), is.character(effectsize), effectsize %in% c("COR","SMD","LOGOR","ZCOR"))
  
  effectsize <- match.arg(effectsize, c("COR","SMD","LOGOR","ZCOR"))
  
  # Depending on the format of the input (correlations, smd, OR) we convert 
  # them to the effect size fisher z transformed correlation for the analysis
  
  # ----------------------------------------------------------------------------
  # Preliminary step to convert different effect sizes into Transformed Fisher Z
  # ----------------------------------------------------------------------------
  
  if (effectsize == "COR"){
    # Check if correlations and sample sizes are in the right format (numeric)
    stopifnot(is.numeric(r_o), is.numeric(n_o), is.numeric(r_r), is.numeric(n_r))
    
    # Convert r to Fisher Z-Transform using the second order Taylor expansion for delta method (preciser than first order expansion)
    es_o <- atanh(r_o) + (r_o / (2 * n_o))
    es_r <- atanh(r_r) + (r_r / (2 * n_r))
    # Convert Std. Err. of r to Fisher Z-Transform using the second order Taylor expansion (preciser than first order expansion)
    se_o <- sqrt( 1 / n_o + (6 - r_o^2) / (2 * n_o^2))
    se_r <- sqrt( 1 / n_r + (6 - r_r^2) / (2 * n_r^2))
    
    inputs <- data.frame(
      es_o = es_o, se_o = se_o, es_r = es_r, se_r = se_r
    )
    
  } else if (effectsize == "SMD"){
    stopifnot(
      is.numeric(d_o), is.numeric(d_r), is.numeric(se_o), is.numeric(se_r), 
      is.numeric(n1_o), is.numeric(n2_o), is.numeric(n1_r), is.numeric(n2_r)
    )
    
    # Convert Cohen's d to correlation, apply Fisher Z-Transform for approx. normality
    # metafor::escalc(measure="COR", )
    es_o <- r_to_zr( d_to_r(d_o, se_o, n1=n1_o , n2=n2_o)[1], d_to_r(d_o, se_o, n1=n1_o , n2=n2_o)[2])[["Fisher Z correlation"]]
    es_r <- r_to_zr( d_to_r(d_r, se_r, n1=n1_r , n2=n2_r)[1], d_to_r(d_r, se_r, n1=n1_r , n2=n2_r)[2])[["Fisher Z correlation"]]
    # Derive standard errors for fisher z Transform using full formula from Stuart and Ord (1994)
    se_o <- r_to_zr( d_to_r(d_o, se_o, n1=n1_o , n2=n2_o)[1], d_to_r(d_o, se_o, n1=n1_o , n2=n2_o)[2])[["Std. Err."]]
    se_r <- r_to_zr( d_to_r(d_r, se_r, n1=n1_r , n2=n2_r)[1], d_to_r(d_r, se_r, n1=n1_r , n2=n2_r)[2])[["Std. Err."]]
    
    inputs <- data.frame(
      es_o = es_o,  se_o = se_o, es_r = es_r, se_r = se_r,
      n1_o = n1_o,  n2_o = n2_o, n1_r = n1_r, n2_r = n2_r
    )
    
  } else if (effectsize == "LOGOR"){
    stopifnot(
      is.numeric(logor_o), is.numeric(logor_r), is.numeric(se_o), is.numeric(se_r),
      is.numeric(n1_o), is.numeric(n2_o), is.numeric(n1_r), is.numeric(n2_r)
    )
    
    # Convert log OR and SE(log OR) to Cohen's d
    d_o <- logor_to_d( logor = logor_o, se_logor = se_o)[["Cohen's d"]]
    d_r <- logor_to_d( logor = logor_r, se_logor = se_r)[["Cohen's d"]]
    se_o <- logor_to_d( logor = logor_o, se_logor = se_o)[["Std. Err."]]
    se_r <- logor_to_d( logor = logor_r, se_logor = se_r)[["Std. Err."]]
    
    # Convert Cohen's d to correlation, apply Fisher Z-Transform for approx. normality
    es_o <- r_to_zr( d_to_r(d_o, se_o, n1=n1_o , n2=n2_o)[1], d_to_r(d_o, se_o, n1=n1_o , n2=n2_o)[2])[["Fisher Z correlation"]]
    es_r <- r_to_zr( d_to_r(d_r, se_r, n1=n1_r , n2=n2_r)[1], d_to_r(d_r, se_r, n1=n1_r , n2=n2_r)[2])[["Fisher Z correlation"]]
    # Derive standard errors for fisher z Transform using full formula from Stuart and Ord (1994)
    se_o <- r_to_zr( d_to_r(d_o, se_o, n1=n1_o , n2=n2_o)[1], d_to_r(d_o, se_o,  n1=n1_o , n2=n2_o)[2])[["Std. Err."]]
    se_r <- r_to_zr( d_to_r(d_r, se_r, n1=n1_r , n2=n2_r)[1], d_to_r(d_r, se_r, n1=n1_r , n2=n2_r)[2])[["Std. Err."]]
    
    # Insert converted effect sizes in a data frame for the analysis
    inputs <- data.frame(
      es_o = es_o, se_o = se_o, es_r = es_r, se_r = se_r,
      n1_o = n1_o, n2_o = n2_o, n1_r = n1_r, n2_r = n2_r
    )
    
  } else if (effectsize == "ZCOR"){
    stopifnot(
      is.numeric(z_o), is.numeric(z_r), is.numeric(se_o), is.numeric(se_r)
    )
    # Insert Fisher Z effect sizes in a data frame for the analysis
    inputs <- data.frame(
      es_o = z_o, se_o = se_o, es_r = z_r, se_r = se_r
    )
  }
  
  # Adding study id for the output table
  inputs$study_id <- seq_len(nrow(inputs))
  
  # ---------------------------------------------------------------------------
  # Now start with the success criteria for each metric
  # ---------------------------------------------------------------------------
  
  # Significance criterion (vote counting) ----
  assess_rep_significance <- function(df, alpha = 0.05) {
    
    #should sig test for replication effects here be one-tailed?
    #LR: should inconsistent direction be always included here? I'd say yes.
    p_o <- zr_to_p(zr = df$es_o, se = df$se_o)
    p_r <- zr_to_p(zr = df$es_r, se = df$se_r)
    
    ### Success criterion  
    
    outcome_detailed <- ifelse(p_o >= alpha,
                               "original effect is not significant",
                               ifelse(p_r >= alpha,
                                      "replication effect is not significant",
                                      ifelse(p_r < alpha & sign(df$es_o) == sign(df$es_r),
                                             "replication effect is significant and same direction of original",
                                             ifelse(p_r < alpha & sign(df$es_o) != sign(df$es_r),
                                                    "direction of effect is inconsistent",
                                                    NA_character_))))
    
    outcome <- ifelse(outcome_detailed == "original effect is not significant", "OS not significant",
                      ifelse(outcome_detailed %in% c("replication effect is not significant",
                                                     "direction of effect is inconsistent"), "Replication failure",
                             ifelse(outcome_detailed == "replication effect is significant and same direction of original", "Replication success",
                                    NA_character_)))
    
    data.frame(outcome=outcome, outcome_detailed=outcome_detailed, quantitative_value = p_r)
  }
  
  
  # Fixed Meta-Analysis (aggregating one original and one replication) ----
  assess_agg_significance <- function(df) {
    
    # Calculates the pooled effect estimate and CI lower bound
    effect_comparison <- compare_effectsizes(df$es_o, df$se_o, df$es_r, df$se_r)
    
    ### Calculate p-value for this metric
    p_meta <- effect_comparison$pval
  
    ### Success criterion  
    outcome_detailed <- ifelse(effect_comparison$ci_lb > 0 & effect_comparison$est > 0,
                               "aggregated effect is larger than 0 and in the same direction",
                               ifelse(effect_comparison$ci_lb > 0 & effect_comparison$est < 0,
                                      "aggregated effect is larger than 0 but in the opposite direction",
                                      ifelse(effect_comparison$ci_lb <= 0,
                                             "aggregated effect is not larger than 0",
                                             NA_character_)))
    
    outcome <- ifelse(outcome_detailed == "aggregated effect is larger than 0 and in the same direction",
                      "Replication success",
                      ifelse(outcome_detailed %in% c("aggregated effect is larger than 0 but in the opposite direction",
                                                     "aggregated effect is not larger than 0"),
                             "Replication failure", NA_character_))
    
    data.frame(outcome=outcome, outcome_detailed=outcome_detailed, quantitative_value = p_meta)
    
  }
  
  # Original effect in Replication Confidence interval metric ----
  assess_or_consistency_with_rep <- function(df) {
    
    ### Calculate p-value for this metric
    
    p_cir <- pnorm(abs(df$es_o - df$es_r) / df$se_r, mean = 0, sd = 1, lower.tail = FALSE)
    
    ### Calculate replication confidence interval
    
    cis <- compute_ci(df$es_r, df$se_r)
    
    ### Success criterion  
    
    outcome_detailed <- ifelse(df$es_o > cis$lower & df$es_o < cis$upper,
                               "Original effect size is within replication's CI",
                               ifelse(df$es_o <= cis$lower | df$es_o >= cis$upper,
                                      "Original effect size is not within replication's CI",
                                      NA_character_))
    
    outcome <- ifelse(outcome_detailed == "Original effect size is within replication's CI",
                      "Replication success",
                      ifelse(outcome_detailed == "Original effect size is not within replication's CI",
                             "Replication failure", NA_character_))
    
    data.frame(outcome=outcome, outcome_detailed=outcome_detailed, quantitative_value = p_cir)
  }
  
  # Replication in Original Confidence interval metric ----
  assess_rep_consistency_with_or <- function(df) {
    
    ### Calculate p-value for this metric
    
    p_cio <- pnorm(abs(df$es_r - df$es_o) / df$se_o, mean = 0, sd = 1, lower.tail = FALSE)
    
    ### Calculate original confidence interval
    cis <- compute_ci(df$es_o, df$se_o)
    
    ### Success criterion  
    
    outcome_detailed <- ifelse(df$es_r > cis$lower & df$es_r < cis$upper,
                               "Replication effect size is within original's CI",
                               ifelse(df$es_r <= cis$lower | df$es_r >= cis$upper,
                                      "Replication effect size is not within original's CI",
                                      NA_character_))
    
    outcome <- ifelse(outcome_detailed == "Replication effect size is within original's CI",
                      "Replication success",
                      ifelse(outcome_detailed == "Replication effect size is not within original's CI",
                             "Replication failure", NA_character_))
    
    data.frame(outcome=outcome, outcome_detailed=outcome_detailed, quantitative_value = p_cio)
  }
  
  # Replication in Original Prediction Interval ----
  assess_rep_consistency_with_pi <- function(df) {
    
    ### Calculate p-value for this metric
    p_PI <- pnorm(abs(df$es_o - df$es_r) / sqrt(df$se_o^2 + df$se_r^2), mean = 0, sd = 1, lower.tail = FALSE)
    
    ### Calculate prediction interval
    pred_intervals <- compute_prediction_interval(df$es_o, df$se_o, df$se_r)
    
    ### Success criterion  
    
    outcome_detailed <- ifelse(df$es_r > pred_intervals$lower & df$es_r < pred_intervals$upper,
                               "Replication effect size is within the prediction interval of the original",
                               ifelse(df$es_r <= pred_intervals$lower | df$es_r >= pred_intervals$upper,
                                      "Replication effect size is not within the prediction interval of the original",
                                      NA_character_))
    
    outcome <- ifelse(outcome_detailed == "Replication effect size is within the prediction interval of the original",
                      "Replication success",
                      ifelse(outcome_detailed == "Replication effect size is not within the prediction interval of the original",
                             "Replication failure", NA_character_))
    
    data.frame(outcome=outcome, outcome_detailed=outcome_detailed, quantitative_value = p_PI)
  }
  
  # Small telescopes Metric ----
  assess_small_telescopes <- function(df, power=0.33, alpha=0.05) {
    if (effectsize %in% c("SMD","LOGOR")){
      es_with_33power <- mapply(find_z_for_power, n1 = df$n1_o, n2 = df$n2_o,
                                MoreArgs = list(power = power))
    } else {
      es_with_33power <- mapply(function(n) {
        power.t.test(power = power, n = n, sig.level = alpha)$delta
      }, n_o) 
    }
    
    p_telescopes <- pnorm(abs(es_with_33power - df$es_r) / df$se_r, mean = 0, sd = 1, lower.tail = FALSE)
    
    ### Success criterion  
    
    outcome_detailed <- ifelse(df$es_r >= es_with_33power,
                               "original study had >= 33% power to detect replication effect",
                               ifelse(df$es_r < es_with_33power,
                                      "original study had < 33% power to detect replication effect",
                                      NA_character_))
    
    outcome <- ifelse(outcome_detailed == "original study had >= 33% power to detect replication effect",
                      "Replication success",
                      ifelse(outcome_detailed == "original study had < 33% power to detect replication effect",
                             "Replication failure", NA_character_))
    
    data.frame(outcome=outcome, outcome_detailed=outcome_detailed, quantitative_value = p_telescopes)
  }
  
  # JZS Bayes Factor ----
  assess_BFJZS <- function(df, gamma = 3){
    
    # Computes integrand of the likelihood and the Cauchy prior (JZS prior)
    integrate_jzs <- function(theta, es, se){
      stats::dnorm(x = es, mean = theta, sd = se) *
        stats::dcauchy(theta, location = 0, scale = 1)
    }
    
    # Function for computing vectorized integral. Divide by marginal likelihood under the null to obtain BF
    compute_bfjzs <- function(es, se) {
      log(stats::integrate(integrate_jzs, lower=-Inf, upper=Inf, es=es, se=se)$value /
            stats::dnorm(x=es, mean=0, sd=se))
    }
    
    BFJZS_ori <- mapply(compute_bfjzs, es = df$es_o, se = df$se_o)
    BFJZS_rep <- mapply(compute_bfjzs, es = df$es_r, se = df$se_r)
    
    ### Success criterion  
    
    outcome_detailed <- ifelse(BFJZS_ori < log(gamma) & BFJZS_ori > log(1/gamma),
                               "Original null effect",
                               ifelse(BFJZS_ori > log(gamma) & BFJZS_rep > log(gamma),
                                      "Evidence for replication success",
                                      ifelse(BFJZS_ori > log(gamma) & BFJZS_rep < log(1/gamma),
                                             "Evidence for replication failure",
                                             NA_character_)))
    
    outcome <- ifelse(BFJZS_ori < log(gamma) & BFJZS_ori > log(1/gamma),
                      "Inconclusive due to original null effect",
                      ifelse(BFJZS_ori > log(gamma) & BFJZS_rep > log(gamma),
                             "Replication success",
                             ifelse(BFJZS_ori > log(gamma) & BFJZS_rep < log(1/gamma),
                                    "Replication failure",
                                    NA_character_)))
    
    data.frame(outcome=outcome, outcome_detailed=outcome_detailed, quantitative_value = exp(BFJZS_rep))
  }
  
  # Equality of effect Bayes Factor ----
  assess_BFe <- function(df, gamma = 1, tau=0){
    ## compute sceptical BF and invert it, so that values >1 indicate evidence consistent with original study (replication success)
    BFe <- BayesRep::BFe(to=df$es_o, so=df$se_o, tr=df$es_r, sr=df$se_r, tau=tau)
    
    ### Success criterion  
    
    outcome_detailed <- ifelse(BFe > gamma,
                               "Evidence for replication success",
                               ifelse(BFe <= gamma,
                                      "Evidence for replication failure",
                                      NA_character_))
    
    outcome <- ifelse(BFe > gamma,   "Replication success",
                      ifelse(BFe <= 1/gamma, "Replication failure",
                             "Inconclusive"))
    
    data.frame(outcome=outcome, outcome_detailed=outcome_detailed, quantitative_value = BFe)
  }
  
  # Replication Bayes Factor ----
  
  assess_BFr <- function(df, gamma = 3){
    ## compute sceptical BF and invert it, so that values >1 indicate evidence consistent with original study (replication success)
    BFr <- BayesRep::BFr(to=df$es_o, so=df$se_o, tr=df$es_r, sr=df$se_r)
    
    ### Success criterion  
    
    outcome_detailed <- ifelse(BFr > gamma,       "Evidence for replication success",
                               ifelse(BFr >= 1,           "Weak evidence for replication success",
                                      ifelse(BFr >= 1/gamma,     "Weak evidence for replication failure",
                                             ifelse(BFr <= 1/gamma,     "Evidence for replication failure",
                                                    NA_character_))))
    
    outcome <- ifelse(BFr > gamma,    "Replication success",
                      ifelse(BFr <= 1/gamma, "Replication failure",
                             "Inconclusive"))
    
    data.frame(outcome=outcome, outcome_detailed=outcome_detailed, quantitative_value = BFr)
  }
  
  # Sceptical Bayes Factor ----
  
  assess_BFs <- function(df, gamma = 3){
    ## compute sceptical BF and invert it, so that values >1 indicate evidence consistent with original study (replication success)
    BFs <- 1 / BayesRep::BFs(to=df$es_o, so=df$se_o, tr=df$es_r, sr=df$se_r)
    
    ### Success criterion  
    
    outcome_detailed <- ifelse(is.nan(BFs),
                               "Results of replication is too inconclusive to establish replication success at any level",
                        ifelse(BFs > gamma,    "Evidence for replication success",
                        ifelse(BFs >= 1,        "Weak evidence for replication success",
                        ifelse(BFs >= 1/gamma,  "Weak evidence for replication failure",
                        ifelse(BFs < 1/gamma,   "Evidence for replication failure",
                               NA_character_)))))
    
    outcome <- ifelse(is.nan(BFs),    "Inconclusive",
               ifelse(BFs > gamma,    "Replication success",
               ifelse(BFs <= 1/gamma, "Replication failure",
                                      "Inconclusive")))
    
    data.frame(outcome=outcome, outcome_detailed=outcome_detailed, quantitative_value = BFs)
  }
  
  # Sceptical p-value ----
  
  assess_pScept <- function(df, alpha=0.05){
    ## compute controlled sceptical p-value, so that values < alpha indicate evidence consistent with original study (replication success)
    pScept <- ReplicationSuccess::pSceptical(
      zo = df$es_o / df$se_o,
      zr = df$es_r / df$se_r,
      c  = df$se_o^2 / df$se_r^2,
      type = "controlled",
      alternative = "two.sided"
    )
    
    ### Success criterion  
    
    outcome_detailed <- ifelse(pScept <= alpha,
                               "Evidence for replication success",
                               ifelse(pScept > alpha,
                                      "No evidence for replication success",
                                      NA_character_))
    
    outcome <- ifelse(pScept <= alpha, "Replication success",
                      ifelse(pScept > alpha,  "Replication failure",
                             NA_character_))
    
    data.frame(outcome=outcome, outcome_detailed=outcome_detailed, quantitative_value = pScept)
  }
  
  # Snapshot hybrid
  
  assess_snaphybrid <- function(df, threshold=0.75){
    
    effect_labels <- c("null effect", "small effect", "medium effect", "large effect")
    
    snap_results <- mapply(function(es_o, es_r, se_o, se_r) {
      vals  <- hybrid_snap(es = c(es_o, es_r), se = c(se_o, se_r))
      probs <- unlist(vals)
      idx   <- which.max(probs)
      c(pSnap = unname(probs[idx]), snap_idx = unname(idx))
    }, df$es_o, df$es_r, df$se_o, df$se_r)
    
    # snap_results is a 2 x n matrix, rows = pSnap/snap_idx, cols = studies
    pSnap    <- snap_results[1, ]
    snap_idx <- as.integer(snap_results[2, ])
    
    ### Success criterion  
    
    outcome_detailed <- ifelse(pSnap >= threshold,
                               paste("Evidence for", effect_labels[snap_idx]),
                               "Inconclusive evidence")
    
    outcome <- ifelse(pSnap >= threshold, "Replication success", "Replication failure")
    
    data.frame(outcome=outcome, outcome_detailed=outcome_detailed, quantitative_value = pSnap)
  }
  
  # Q-test ----
  
  assess_Qtest <- function(df, alpha = 0.05){
    # Compute Q-statistic (related to Prediction Interval metric)
    pQtest <- ReplicationSuccess::Qtest(thetao = df$es_o, thetar = df$es_r, 
                                       seo = df$se_o, ser = df$se_r)
    
    ### Success criterion  
    
    outcome_detailed <- ifelse(pQtest <= alpha/2,
                               "Evidence for replication success",
                               ifelse(pQtest > alpha/2,
                                      "No evidence for replication success",
                                      NA_character_))
    
    outcome <- ifelse(pQtest <= alpha/2, "Replication success",
                      ifelse(pQtest > alpha/2,  "Replication failure",
                             NA_character_))
    
    data.frame(outcome=outcome, outcome_detailed=outcome_detailed, quantitative_value = pQtest)
  }
    
  # Likelihood ratio approach
  
  assess_Likratio <- function(df){
    
  }
  
  # #######################################################################----
  # Now create a data frame with the results
  # #######################################################################----
  
  # Switching the names of the functions to metric argument
  metric_fns <- list(
    significance_r   = assess_rep_significance,
    significance_agg = assess_agg_significance,
    consistency_cir  = assess_or_consistency_with_rep,
    consistency_cio  = assess_rep_consistency_with_or,
    consistency_pi   = assess_rep_consistency_with_pi,
    small_telescopes = assess_small_telescopes,
    BF_JZS           = assess_BFJZS,
    BFe              = assess_BFe,
    BFr              = assess_BFr,
    BFs              = assess_BFs,
    pScept           = assess_pScept,
    snaphybrid       = assess_snaphybrid,
    Qtest            = assess_Qtest
  )
  
  # Validate that the inserted metrics are correct
  unknown <- setdiff(metric, names(metric_fns))
  if (length(unknown) > 0) {
    stop("Metric(s) not recognised: ", paste(unknown, collapse=", "),
         "\nAvailable metrics: ", paste(names(metric_fns), collapse=", "))
  }
  
  # Run each metric and collect results
  results <- lapply(metric, function(m) {
    res <- metric_fns[[m]](inputs)
    res$study_id <- inputs$study_id
    res$metric   <- m
    res$es_o     <- inputs$es_o
    res$se_o     <- inputs$se_o
    res$es_r     <- inputs$es_r
    res$se_r     <- inputs$se_r
    # Reorder so identifiers come first
    res[, c("study_id","metric","quantitative_value","es_o","se_o","es_r","se_r",
            "outcome","outcome_detailed")]
  })
  
  output <- do.call(rbind, results)
  output <- output[order(output$study_id, output$metric), ]
  rownames(output) <- NULL
  
  return(output)
}

