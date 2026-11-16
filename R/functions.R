###############################################################################
#############################  ReproMeta  #####################################
###############################################################################
library(ReplicationSuccess)
library(BayesRep)
library(metafor)

# Main idea: Meta-function taking different effect sizes as input.
# Correlation, cohen's d, log OR will be converted to transformed Fisher z.

# %||% is only in base R from 4.4.0; define a fallback for older versions
if (!exists("%||%")) `%||%` <- function(x, y) if (is.null(x)) y else x

# Function 1: Effect size conversion ----------------------------------------

#' prepare_effect_sizes: Effect size conversion
#'
#' Converts effect-size estimates from an original study and its replication
#' to a common format used by the replication-assessment functions. Supported
#' input formats are correlations, standardized mean differences, log odds
#' ratios, and precomputed Fisher-\eqn{z} correlations.
#'
#' @param r_o,r_r Numeric vectors containing the correlations from the original
#'   and replication studies. Required when `effectsize = "COR"`.
#' @param n_o,n_r Numeric vectors containing the sample sizes of the original
#'   and replication studies. Required for `"COR"` and `"ZCOR"` effect sizes.
#' @param effectsize Character vector specifying the effect-size format.
#'   Available options are:
#'   \describe{
#'     \item{`"COR"`}{Correlation coefficients.}
#'     \item{`"SMD"`}{Standardized mean differences.}
#'     \item{`"LOGOR"`}{Log odds ratios.}
#'     \item{`"ZCOR"`}{Precomputed Fisher-\eqn{z} correlations.}
#'   }
#'   A single value is recycled across all study pairs. Alternatively, one
#'   value may be supplied for each pair.
#' @param scale Character string specifying the analysis scale for correlation
#'   inputs. Must be either `"fisherz"` (the default) or `"cor"`.
#'   When `"fisherz"` is selected, correlations and any supplied confidence
#'   interval bounds are transformed using the inverse hyperbolic tangent.
#' @param d_o,d_r Numeric vectors containing standardized mean differences from
#'   the original and replication studies. Required when
#'   `effectsize = "SMD"`.
#' @param se_o,se_r Numeric vectors containing the standard errors of the
#'   supplied effect sizes. For `"SMD"`, these are standard errors of the
#'   standardized mean differences; for `"LOGOR"`, they are standard errors
#'   of the log odds ratios; and for `"ZCOR"`, they are standard errors of the
#'   Fisher-\eqn{z} correlations. For `"COR"`, standard errors are calculated
#'   internally.
#' @param logor_o,logor_r Numeric vectors containing log odds ratios from the
#'   original and replication studies. Required when
#'   `effectsize = "LOGOR"`.
#' @param n1_o,n2_o Numeric vectors containing the group-specific sample sizes
#'   in the original study. Required for `"SMD"` and `"LOGOR"`.
#' @param n1_r,n2_r Numeric vectors containing the group-specific sample sizes
#'   in the replication study. Required for `"SMD"` and `"LOGOR"`.
#' @param z_o,z_r Numeric vectors containing precomputed Fisher-\eqn{z}
#'   correlations from the original and replication studies. Required when
#'   `effectsize = "ZCOR"`.
#' @param orig_conv_r_lb,orig_conv_r_ub Optional numeric vectors containing
#'   lower and upper confidence limits for the original correlations.
#'   The supplied limits must be on the raw correlation scale.
#' @param repli_conv_r_lb,repli_conv_r_ub Optional numeric vectors containing
#'   lower and upper confidence limits for the replication correlations.
#'   The supplied limits must be on the raw correlation scale.
#' @param orig_p_value Optional numeric vector containing reported original-study
#'   p-values. These values may be used by methods such as the sum of p-values
#'   and the sceptical p-value.
#' @param orig_p_value_type Optional character vector describing the reporting
#'   type of each original-study p-value. When supplied, only values coded as
#'   `"exact"` are retained.
#' @param orig_p_value_tails Optional character vector indicating whether each
#'   original-study p-value is one- or two-tailed. Reported one-tailed p-values
#'   are converted to two-tailed p-values by multiplication by two.
#' @param repli_p_value Optional numeric vector containing reported replication
#'   p-values.
#' @param stat_type_o,stat_type_r Optional character vectors identifying the
#'   test-statistic type in the original and replication studies. A value of
#'   `"t"` causes the standard error calculation to use degrees of freedom;
#'   other test types use the effective sample size.
#' @param dof_o,dof_r Optional numeric vectors containing the reported degrees
#'   of freedom in the original and replication studies.
#' @param effective_df_o,effective_df_r Optional numeric vectors containing
#'   effective degrees of freedom. When available, these take precedence over
#'   `dof_o` and `dof_r`.
#' @param effective_sample_size_o,effective_sample_size_r Optional numeric
#'   vectors containing effective sample sizes. When available, these take
#'   precedence over `n_o` and `n_r` in the standard-error
#'   calculations.
#' @param pattern_met Optional logical vector indicating whether the
#'   replication result satisfies the qualitative direction or pattern
#'   criterion. When unavailable, metrics requiring this information may use a
#'   sign-agreement approximation.
#'
#' @details
#' For correlation inputs, the function can retain the raw correlation scale
#' or transform the estimates to the Fisher-\eqn{z} scale. If test-statistic
#' metadata are supplied, standard errors are calculated
#' using effective degrees of freedom or effective sample sizes.
#'
#' Standardized mean differences are first converted to correlations and then,
#' when `scale = "fisherz"`, to Fisher-\eqn{z} correlations. Log odds ratios
#' are converted successively to standardized mean differences, correlations,
#' and (for `scale = "fisherz"`) Fisher-\eqn{z} correlations. Inputs supplied
#' as `"ZCOR"` are retained on their supplied scale.
#'
#' Optional reported confidence intervals, p-values, and pattern indicators
#' are attached to the resulting data frame.
#'
#' @return A data frame with one row per original–replication pair. It contains
#'   at least the following columns:
#'   \describe{
#'     \item{`study_id`}{Identifier for the study pair.}
#'     \item{`effectsize`}{Input effect-size format.}
#'     \item{`scale`}{Requested analysis scale.}
#'     \item{`es_o`, `es_r`}{Prepared original and replication effect sizes.}
#'     \item{`se_o`, `se_r`}{Corresponding standard errors.}
#'     \item{`n_o`, `n_r`}{Total original and replication sample sizes.}
#'   }
#'   Group-specific sample sizes and optional confidence intervals, p-values,
#'   or pattern indicators are included when supplied. Columns containing only
#'   missing values are removed.
#'
#' @export
#'
#' @examples
#' # Prepare correlations on the Fisher-z scale
#' dat <- prepare_effect_sizes(
#'   r_o = c(0.30, 0.20),
#'   n_o = c(100, 150),
#'   r_r = c(0.25, 0.05),
#'   n_r = c(120, 180),
#'   effectsize = "COR",
#'   scale = "fisherz"
#' )
#'
#' dat
#'
#' # Keep correlations on their original scale
#' dat_cor <- prepare_effect_sizes(
#'   r_o = 0.30,
#'   n_o = 100,
#'   r_r = 0.25,
#'   n_r = 120,
#'   effectsize = "COR",
#'   scale = "cor"
#' )
prepare_effect_sizes <- function(r_o=NULL, n_o=NULL, r_r=NULL, n_r=NULL,
                                 effectsize = "COR",
                                 scale = c("fisherz", "cor"),
                                 d_o=NULL, d_r=NULL, se_o=NULL, se_r=NULL,
                                 logor_o=NULL, logor_r=NULL,
                                 n1_o=NULL, n2_o=NULL, n1_r=NULL, n2_r=NULL,
                                 z_o=NULL, z_r=NULL,
                                 # Optional arguments: precomputed CI bounds, ON THE CORRELATION SCALE.
                                 # Attached to the output on the analysis scale: raw when scale="cor",
                                 # atanh()-transformed when scale="fisherz". When present,
                                 # the metrics computing CIs (e.g., consistency_cir/cio,
                                 # small_telescopes) use them instead of recomputing a CI.
                                 orig_conv_r_lb=NULL, orig_conv_r_ub=NULL,
                                 repli_conv_r_lb=NULL, repli_conv_r_ub=NULL,
                                 # Optional arguments: reported p-values + coding metadata.
                                 # orig_p_value is kept only when orig_p_value_type =="exact" (if type is supplied),
                                 # doubled when reported one-tailed (if tails is supplied)
                                 # with limits to (0,1] if correlation.
                                 # These provided p-values are used by sum_p and by pScept via the p2z() route.
                                 orig_p_value=NULL, orig_p_value_type=NULL,
                                 orig_p_value_tails=NULL, repli_p_value=NULL,
                                 # Optional: test-statistic metadata for exact, test-type-dependent SE.
                                 # When supplied, se_o/se_r are overwritten: on the cor scale with
                                 # sqrt((1-r^2)/df1) for t-type (else /size); on the
                                 # fisherz scale with 1/sqrt(n_eff - 3), where n_eff is
                                 # the same effective df/size.
                                 stat_type_o=NULL, dof_o=NULL,
                                 effective_df_o=NULL, effective_sample_size_o=NULL,
                                 stat_type_r=NULL, dof_r=NULL,
                                 effective_df_r=NULL, effective_sample_size_r=NULL,
                                 # Optional: qualitative direction/pattern flag.
                                 # Used to gate pScept / BFr on the cor scale
                                 # Without this argument these metrics fall back to
                                 # a sign-match approximation.
                                 pattern_met=NULL,
                                 # Optional: TRUE if the replication effect is in the direction of the
                                 # original claim (e.g., SCORE's repli_effect_direction). Used by the
                                 # metrics that combine p-values on the Fisher-z scale.
                                 same_direction=NULL) {

  scale <- match.arg(scale)
  stopifnot(is.character(effectsize), all(effectsize %in% c("COR","SMD","LOGOR","ZCOR")))

  n_studies <- length(r_o %||% d_o %||% logor_o %||% z_o)
  if (length(effectsize) == 1) effectsize <- rep(effectsize, n_studies)
  stopifnot(length(effectsize) == n_studies)

  # Exact test-type-dependent SE on the raw correlation scale
  # (sqrt((1-r^2)/df1) if stat_type=="t", else sqrt((1-r^2)/size),
  #  preferring effective df / effective sample size, falling back to raw dof/n)
  compute_se_score <- function(r, stat_type, dof, effective_df,
                               sample_size, effective_sample_size) {
    if (is.na(r) || is.na(stat_type)) return(NA_real_)
    df1  <- ifelse(is.na(effective_df), dof, effective_df)
    size <- ifelse(is.na(effective_sample_size), sample_size, effective_sample_size)
    if (stat_type == "t") {
      if (is.na(df1) || df1 <= 0) return(NA_real_)
      return(sqrt((1 - r^2) / df1))
    }
    if (is.na(size) || size <= 0) return(NA_real_)
    sqrt((1 - r^2) / size)
  }

  # Fisher-z SE, in order of availability:
  #  1. degrees of freedom of a t test (effective df if reported, else dof):
  #     df1 = n - 2 - k with k covariates and Var(atanh(r)) = 1 / (n - 3 - k),
  #     so se = 1 / sqrt(df1 - 1).
  #  2. effective sample size: se = 1 / sqrt(n_eff - 3)
  #  3. sample size (the default without additional data, default when additional arguments are empty)
  compute_se_z <- function(r, n, stat_type, dof, effective_df, effective_sample_size) {
    df1 <- ifelse(is.na(effective_df), dof, effective_df)
    if (!is.na(stat_type) && stat_type == "t" && !is.na(df1) && df1 > 1)
      return(list(se = 1 / sqrt(df1 - 1), src = "df"))
    if (!is.na(effective_sample_size) && effective_sample_size > 3)
      return(list(se = 1 / sqrt(effective_sample_size - 3), src = "effective_n"))
    list(se = sqrt(1 / n + (6 - r^2) / (2 * n^2)), src = "n")
  }

  # Cohen's d (and its SE) -> point-biserial r, then Fisher z unless scale = "cor"
  d_to_analysis_scale <- function(d, se_d, n1, n2, scale) {
    r <- d_to_r(d, se_d, n1 = n1, n2 = n2)
    if (scale == "cor") return(c(es = r[["r"]], se = r[["Std. Err."]]))
    z <- r_to_zr(r[["r"]], r[["Std. Err."]])
    c(es = z[["Fisher Z correlation"]], se = z[["Std. Err."]])
  }

  # Checks if the SE-adjustment metadata supplied at all (e.g., if provided p-values are two-sided resp. one-sided they need to be adjusted)
  have_se_meta <- !is.null(stat_type_o) || !is.null(stat_type_r) ||
    !is.null(dof_o) || !is.null(dof_r) ||
    !is.null(effective_df_o) || !is.null(effective_df_r) ||
    !is.null(effective_sample_size_o) || !is.null(effective_sample_size_r)

  # pad any supplied-but-missing metadata vector to full length with NA
  pad <- function(vec) if (is.null(vec)) rep(NA, n_studies) else vec

  convert_one <- function(i) {
    es <- effectsize[i]

    if (es == "COR") {
      if (is.null(r_o) || is.null(n_o) || is.null(r_r) || is.null(n_r))
        stop("effectsize 'COR' requires r_o, n_o, r_r, n_r - missing for study ", i)

      src_o <- src_r <- "n"
      if (scale == "cor") {
        es_o_i <- r_o[i]
        es_r_i <- r_r[i]
        if (have_se_meta) {
          src_o <- src_r <- "test"
          # test-type-dependent SE
          se_o_i <- compute_se_score(r_o[i], pad(stat_type_o)[i], pad(dof_o)[i],
                                     pad(effective_df_o)[i], n_o[i],
                                     pad(effective_sample_size_o)[i])
          se_r_i <- compute_se_score(r_r[i], pad(stat_type_r)[i], pad(dof_r)[i],
                                     pad(effective_df_r)[i], n_r[i],
                                     pad(effective_sample_size_r)[i])
        } else {
          se_o_i <- sqrt((1 - r_o[i]^2) / (n_o[i] - 3))
          se_r_i <- sqrt((1 - r_r[i]^2) / (n_r[i] - 3))
        }
      } else {
        # Exact Fisher-z transform containment checks against atanh()-transformed
        # CI bounds only if es is atanh(r) exactly.
        es_o_i <- atanh(r_o[i])
        es_r_i <- atanh(r_r[i])
        # degrees of freedom -> effective sample size -> sample size; without
        # test metadata this is the sample-size SE, as before
        se_o_l <- compute_se_z(r_o[i], n_o[i], pad(stat_type_o)[i], pad(dof_o)[i],
                               pad(effective_df_o)[i], pad(effective_sample_size_o)[i])
        se_r_l <- compute_se_z(r_r[i], n_r[i], pad(stat_type_r)[i], pad(dof_r)[i],
                               pad(effective_df_r)[i], pad(effective_sample_size_r)[i])
        se_o_i <- se_o_l$se; src_o <- se_o_l$src
        se_r_i <- se_r_l$se; src_r <- se_r_l$src
      }

      data.frame(es_o=es_o_i, se_o=se_o_i, es_r=es_r_i, se_r=se_r_i,
                 n_o=n_o[i], n_r=n_r[i],
                 n1_o=NA_real_, n2_o=NA_real_, n1_r=NA_real_, n2_r=NA_real_,
                 se_src_o=src_o, se_src_r=src_r)

    } else if (es == "SMD") {
      if (is.null(d_o) || is.null(d_r) || is.null(se_o) || is.null(se_r) ||
          is.null(n1_o) || is.null(n2_o) || is.null(n1_r) || is.null(n2_r))
        stop("effectsize 'SMD' requires d_o, d_r, se_o, se_r, n1_o, n2_o, n1_r, n2_r - missing for study ", i)

      conv_o <- d_to_analysis_scale(d_o[i], se_o[i], n1_o[i], n2_o[i], scale)
      conv_r <- d_to_analysis_scale(d_r[i], se_r[i], n1_r[i], n2_r[i], scale)
      es_o_i <- conv_o[["es"]]; se_o_i <- conv_o[["se"]]
      es_r_i <- conv_r[["es"]]; se_r_i <- conv_r[["se"]]

      data.frame(es_o=es_o_i, se_o=se_o_i, es_r=es_r_i, se_r=se_r_i,
                 n_o=n1_o[i]+n2_o[i], n_r=n1_r[i]+n2_r[i],
                 n1_o=n1_o[i], n2_o=n2_o[i], n1_r=n1_r[i], n2_r=n2_r[i],
                 se_src_o="supplied", se_src_r="supplied")

    } else if (es == "LOGOR") {
      if (is.null(logor_o) || is.null(logor_r) || is.null(se_o) || is.null(se_r) ||
          is.null(n1_o) || is.null(n2_o) || is.null(n1_r) || is.null(n2_r))
        stop("effectsize 'LOGOR' requires logor_o, logor_r, se_o, se_r, n1_o, n2_o, n1_r, n2_r - missing for study ", i)

      d_o_i  <- logor_to_d(logor=logor_o[i], se_logor=se_o[i])[["Cohen's d"]]
      d_r_i  <- logor_to_d(logor=logor_r[i], se_logor=se_r[i])[["Cohen's d"]]
      se_o_i <- logor_to_d(logor=logor_o[i], se_logor=se_o[i])[["Std. Err."]]
      se_r_i <- logor_to_d(logor=logor_r[i], se_logor=se_r[i])[["Std. Err."]]

      conv_o <- d_to_analysis_scale(d_o_i, se_o_i, n1_o[i], n2_o[i], scale)
      conv_r <- d_to_analysis_scale(d_r_i, se_r_i, n1_r[i], n2_r[i], scale)
      es_o_i <- conv_o[["es"]]; se_o_i <- conv_o[["se"]]
      es_r_i <- conv_r[["es"]]; se_r_i <- conv_r[["se"]]

      data.frame(es_o=es_o_i, se_o=se_o_i, es_r=es_r_i, se_r=se_r_i,
                 n_o=n1_o[i]+n2_o[i], n_r=n1_r[i]+n2_r[i],
                 n1_o=n1_o[i], n2_o=n2_o[i], n1_r=n1_r[i], n2_r=n2_r[i],
                 se_src_o="supplied", se_src_r="supplied")

    } else if (es == "ZCOR") {
      if (is.null(z_o) || is.null(z_r) || is.null(se_o) || is.null(se_r) || is.null(n_o) || is.null(n_r))
        stop("effectsize 'ZCOR' requires z_o, z_r, se_o, se_r, n_o, n_r - missing for study ", i)

      data.frame(es_o = z_o[i], se_o = se_o[i], es_r = z_r[i], se_r = se_r[i],
                 n_o = n_o[i], n_r = n_r[i],
                 n1_o = NA_real_, n2_o = NA_real_, n1_r = NA_real_, n2_r = NA_real_,
                 se_src_o = "supplied", se_src_r = "supplied")
    }
  }

  df <- do.call(rbind, lapply(seq_len(n_studies), convert_one))
  df <- cbind(study_id=seq_len(n_studies), effectsize=effectsize, scale=scale, df)
  rownames(df) <- NULL
  df <- df[, colSums(is.na(df)) < nrow(df)]

  # Attach precomputed CI bounds on the analysis scale
  # Input bounds are on the correlation scale
  # scale="cor"    -> attached raw
  # scale="fisherz"-> atanh()-transformed (bounds clamped just inside (-1,1) so
  #                   atanh stays finite), under es_lb_*/es_ub_* names.
  # CI-based assessments against the transformed bounds are identical to raw-r decisions.
  # A side that is supplied without its partner is attached with the partner
  # as NA, so that the metrics use the reported side and compute the other.
  if (!is.null(orig_conv_r_lb) || !is.null(orig_conv_r_ub) ||
      !is.null(repli_conv_r_lb) || !is.null(repli_conv_r_ub)) {
    z_bound <- function(b) atanh(pmax(pmin(b, 1 - 1e-12), -1 + 1e-12))
    to_scale <- if (scale == "cor") identity else z_bound
    side <- function(x, partner) {
      if (!is.null(x)) x else if (!is.null(partner)) rep(NA_real_, n_studies) else NULL
    }
    bounds <- list(side(orig_conv_r_lb, orig_conv_r_ub), side(orig_conv_r_ub, orig_conv_r_lb),
                   side(repli_conv_r_lb, repli_conv_r_ub), side(repli_conv_r_ub, repli_conv_r_lb))
    names(bounds) <- if (scale == "cor") {
      c("orig_conv_r_lb", "orig_conv_r_ub", "repli_conv_r_lb", "repli_conv_r_ub")
    } else {
      c("es_lb_o", "es_ub_o", "es_lb_r", "es_ub_r")
    }
    for (nm in names(bounds)) if (!is.null(bounds[[nm]])) df[[nm]] <- to_scale(bounds[[nm]])
  }

  #  Attach standardized p-values (scale-free)
  #  one-tailed original p-values are doubled to two-tailed
  if (!is.null(orig_p_value)) {
    p <- orig_p_value
    if (!is.null(orig_p_value_type)) {
      p <- ifelse(is.na(orig_p_value_type) | orig_p_value_type != "exact",
                  NA_real_, p)
    }
    if (!is.null(orig_p_value_tails)) {
      p <- ifelse(is.na(orig_p_value_tails) | grepl("two", orig_p_value_tails),
                  p, p * 2)
    }
    df$orig_p_value <- pmin(pmax(p, 1e-10), 1)
  }
  if (!is.null(repli_p_value)) {
    df$repli_p_value <- pmin(pmax(repli_p_value, 1e-10), 1)
  }

  # Qualitative direction/pattern flag
  if (!is.null(pattern_met)) {
    df$repli_pattern_criteria_met <- pattern_met
  }

  # Direction of the replication relative to the original claim
  if (!is.null(same_direction)) {
    df$same_direction <- as.logical(same_direction)
  }

  return(df)
}


# Function 2: Assess replication success criteria ----

#' assess_replication_outcome: Assess replication success criteria
#'
#' Applies one or more replication-success metrics to pairs of original and
#' replication effect sizes prepared with [prepare_effect_sizes()].
#'
#' @param df A data frame containing the prepared original and replication
#'   results. It must contain the columns `study_id`, `es_o`, `se_o`, `es_r`,
#'   and `se_r`. Most metrics additionally require `scale`, `effectsize`,
#'   `n_o`, and `n_r`. Optional reported p-values, confidence intervals, and
#'   qualitative pattern indicators are used when available.
#' @param metric Character vector specifying the replication-success metrics
#'   to calculate. Available metrics are:
#'   \describe{
#'     \item{`"significance_r"`}{Two-trials rule: both studies must provide
#'       statistically significant results in the same direction.}
#'     \item{`"significance_agg"`}{Significance of the effect obtained by
#'       combining the original and replication evidence.}
#'     \item{`"sum_p"`}{Sum-of-p-values criterion with agreement in the
#'       qualitative effect pattern.}
#'     \item{`"consistency_cir"`}{Whether the original effect lies within the
#'       replication confidence interval.}
#'     \item{`"consistency_cio"`}{Whether the replication effect lies within
#'       the original confidence interval.}
#'     \item{`"consistency_pi"`}{Whether the replication effect lies within
#'       the prediction interval derived from the original study.}
#'     \item{`"small_telescopes"`}{Small-telescopes criterion based on the
#'       effect that the original study had 33 percent power to detect.}
#'     \item{`"BF_JZS"`}{Independent Jeffreys-Zellner-Siow Bayes-factor
#'       criterion: both Bayes factors above the threshold and the two
#'       estimates in the same direction.}
#'     \item{`"BF_meta"`}{Fixed-effect meta-analytic Bayes factor: JZS Bayes
#'       factor for an effect common to the original and replication studies,
#'       computed from their inverse-variance pooled estimate; success also
#'       requires the pooled estimate to have the direction of the original.}
#'     \item{`"BFe"`}{Equality-of-effect-sizes Bayes factor.}
#'     \item{`"BFr"`}{Replication Bayes factor.}
#'     \item{`"BFs"`}{Sceptical Bayes factor.}
#'     \item{`"pScept"`}{Sceptical p-value.}
#'     \item{`"snaphybrid"`}{Snapshot hybrid based on posterior probabilities
#'       for the null, small, medium, and large effect models.}
#'     \item{`"equivalence"`}{Equivalence test for the difference between the
#'       original and replication effects.}
#'     \item{`"correspondence"`}{Correspondence test distinguishing
#'       equivalence, trivial difference, meaningful difference, and
#'       indeterminacy.}
#'     \item{`"Likratio"`}{Likelihood-ratio comparison of a theoretically
#'       interesting effect against the null effect.}
#'     \item{`"minimum_effect"`}{Minimum-effect test for a meaningful
#'       difference between the original and replication effects.}
#'   }
#' @param calibrated Logical scalar indicating whether calibrated success
#'   thresholds should be used. The default is `TRUE`. Calibration targets an
#'   overall T1E rate of \eqn{\alpha^2} (by default
#'   \eqn{0.025^2 = 0.000625}) for metrics for which a calibrated version is
#'   implemented. The thresholds are computed from `alpha` with
#'   [calibrate_threshold()]:
#'   \itemize{
#'     \item `"significance_agg"` uses an unweighted Stouffer combination with
#'       threshold \eqn{\sqrt{2} z_{1-\alpha^2}};
#'     \item `"BF_JZS"` uses the Bayes-factor threshold
#'       `calibrate_threshold("BF_JZS", alpha)` (1.60 for alpha = 0.025), for
#'       which the overall T1E rate is at most \eqn{\alpha^2} (for a
#'       pre-specified direction) for every standard error;
#'     \item `"BF_meta"` uses the Bayes-factor threshold
#'       `calibrate_threshold("BF_meta", alpha)` (23.7 for alpha = 0.025),
#'       for which the overall T1E rate is at most \eqn{\alpha^2} for every
#'       standard error of the pooled estimate;
#'     \item `"BFs"` uses the Bayes-factor threshold
#'       `calibrate_threshold("BFs", alpha, c)` computed from the variance
#'       ratio \eqn{c = v_o / v_r} of each pair (2.81 for
#'       alpha = 0.025 and c = 1), which gives an overall T1E rate of exactly
#'       \eqn{\alpha^2} for every pair; and
#'     \item `"pScept"` uses the controlled rather than the golden
#'       recalibration at level `alpha`.
#'   }
#'   The two-trials rule and the sum-of-p-values criterion need no
#'   recalibration: at one-sided level `alpha` their overall T1E rate is
#'   already `alpha^2` (for a pre-specified direction), so they use `alpha`
#'   whatever the value of `calibrated`. Family 2 metrics and the snapshot
#'   hybrid are not calibrated and ignore `alpha`: they keep their own
#'   thresholds (95% intervals, TOST at 0.05, gamma = 1, lr = 8, posterior
#'   probability 0.75), which can be changed with `metric_args`.
#' @param alpha One-sided significance level of the Family 1 metrics
#'   (`"significance_r"`, `"significance_agg"`, `"sum_p"`, `"pScept"`, and,
#'   when `calibrated = TRUE`, `"BF_JZS"`, `"BF_meta"` and `"BFs"`). The overall T1E rate
#'   targeted by the calibrated versions is `alpha^2`. Default 0.025.
#'   Family 2 metrics keep their own conventional levels.
#' @param metric_args Optional named list with one element per metric, each a
#'   named list of arguments that replace that metric's defaults, e.g.
#'   `list(BFe = list(gamma = 3), consistency_pi = list(alpha = 0.05))`.
#'   `NULL` (default) uses the defaults everywhere. The available arguments
#'   are:
#'   \describe{
#'     \item{`consistency_cir`, `consistency_cio`, `consistency_pi`}{`alpha`:
#'       one-sided level, the interval is the (1 - 2 `alpha`) interval
#'       (default 0.025, i.e. 95%). Reported confidence limits (correlation
#'       scale, or `es_lb_*`/`es_ub_*` on the Fisher-z scale) are rescaled
#'       on the Fisher-z scale around the point estimate when `alpha`
#'       differs from 0.025.}
#'     \item{`small_telescopes`}{`alpha`: level of the test against the
#'       benchmark effect (default 0.025; reported limits are rescaled as
#'       above); `alpha_benchmark`: level used to define the effect with
#'       33 percent power in the original study (default `alpha`); `power`
#'       (default 1/3).}
#'     \item{`BF_JZS`, `BF_meta`, `BFs`}{`gamma`: evidence threshold
#'       (default 1 when `calibrated = FALSE`). A value given here replaces
#'       the calibrated threshold.}
#'     \item{`BFr`}{`gamma` (default 1).}
#'     \item{`BFe`}{`gamma` (default 1), `tau` (prior standard deviation,
#'       default 0.2).}
#'     \item{`snaphybrid`}{`threshold`: posterior probability the winning
#'       model must exceed (default 0.75); prior probabilities `pr_0`,
#'       `pr_s`, `pr_m` (default 0.25 each).}
#'     \item{`Likratio`}{`lr`: evidence threshold (default 8); `lr_dtie`:
#'       likelihood ratio that defines the theoretically interesting effect
#'       (default `lr`).}
#'     \item{`equivalence`, `correspondence`}{`h` (default 0.2, correlation
#'       units), `alpha` (default 0.05).}
#'     \item{`minimum_effect`}{`h` (default `NULL`, i.e.
#'       (1 - `shrinkage`) times the original effect), `shrinkage`
#'       (default 0.25), `alpha` (default 0.025).}
#'     \item{`significance_r`, `significance_agg`, `sum_p`, `pScept`}{`alpha`
#'       (default: the `alpha` argument above, which is the recommended way
#'       to change it); `sum_p` also takes `threshold`.}
#'   }
#'
#' @details
#' Each selected metric is applied separately to every original-replication
#' pair. The function returns both a binary or categorical replication outcome
#' and a more detailed interpretation.
#'
#' The interpretation of the quantitative value depends on the selected
#' metric. It may be a p-value, Bayes factor, posterior model probability,
#' likelihood ratio, effect-size threshold, or a missing value for metrics
#' defined only through interval containment.
#'
#' When specific information is available in `df`, the function uses:
#' \itemize{
#'   \item reported confidence limits for the interval-containment and
#'     small-telescopes metrics;
#'   \item exact reported p-values for the sum-of-p-values and sceptical
#'     p-value metrics;
#'   \item `repli_pattern_criteria_met` to assess qualitative agreement; and
#'   \item `snaphybrid_base_eligible` to reproduce the analysis sample used
#'     for the snapshot-hybrid method.
#' }
#'
#' If a required value is missing, invalid, or the metric cannot be calculated,
#' the corresponding output is generally returned as `NA`.
#'
#' @return A named list with two data frames:
#'   \describe{
#'     \item{`data`}{The input effect-size information, containing `study_id`,
#'       `es_o`, `se_o`, `es_r`, and `se_r`.}
#'     \item{`results`}{One row per study pair in wide format. For every
#'       selected metric `m`, the data frame contains:
#'       \describe{
#'         \item{`m_type`}{Label describing the type of quantitative value.}
#'         \item{`m_value`}{Metric-specific quantitative value.}
#'         \item{`m_outcome`}{Overall replication classification.}
#'         \item{`m_outcome_detailed`}{Detailed interpretation of the result.}
#'       }
#'     }
#'   }
#'
#' @seealso [prepare_effect_sizes()]
#'
#' @export
#'
#' @examples
#' dat <- prepare_effect_sizes(
#'   r_o = c(0.30, 0.20),
#'   n_o = c(100, 150),
#'   r_r = c(0.25, 0.05),
#'   n_r = c(120, 180),
#'   effectsize = "COR",
#'   scale = "fisherz"
#' )
#'
#' # Apply one replication metric
#' assess_replication_outcome(
#'   dat,
#'   metric = "significance_r"
#' )
#'
#' # Apply several metrics
#' assess_replication_outcome(
#'   dat,
#'   metric = c(
#'     "significance_r",
#'     "significance_agg",
#'     "consistency_pi",
#'     "BF_JZS"
#'   ),
#'   calibrated = TRUE
#' )
#'
#' # Change thresholds of individual metrics
#' assess_replication_outcome(
#'   dat,
#'   metric = c("consistency_pi", "BFe", "Likratio"),
#'   metric_args = list(
#'     consistency_pi = list(alpha = 0.05),   # 90% prediction interval
#'     BFe            = list(gamma = 3),
#'     Likratio       = list(lr = 10)
#'   )
#' )
assess_replication_outcome <- function(df, metric, calibrated = TRUE, alpha = 0.025,
                                       metric_args = NULL) {

  # `calibrated = TRUE` switches conventional thresholds that do NOT control
  # the overall T1E rate at alpha^2 (default 0.025^2)
  # over to calibrated versions, with thresholds computed from `alpha` by
  # calibrate_threshold():
  #   BF_JZS  : gamma 1 -> calibrate_threshold("BF_JZS", alpha)          (1.60 for alpha = 0.025; overall T1E <= 0.025^2 for every SE)
  #   BFs     : gamma 1 -> calibrate_threshold("BFs", alpha, c) per pair (2.36 for alpha = 0.025 and c = 1; overall T1E = 0.025^2 for every c)
  #   BF_meta : gamma 1 -> calibrate_threshold("BF_meta", alpha)         (23.7 for alpha = 0.025; overall T1E <= 0.025^2 for every SE)
  #   pScept  : golden -> 'controlled' recalibration at level alpha^2
  #   meta    : fixed-effect pooling -> unweighted Stouffer at sqrt(2)*z_{1-alpha^2}
  # Metrics already calibrated at the alpha^2 (two-trials rule,
  # sum of p-values) use `alpha` directly.
  # Since some thresholds (such as gamma, lr) are defined inside metric-specific functions,
  # setting alpha=0.025 in assess_replication_outcome won't pass it down.
  # `metric_args` passes metric-specific arguments (thresholds such as gamma,
  # lr, the snapshot-hybrid threshold, or the alpha of Family 2 metrics) on to
  # the metric functions below. If an explicit gamma for BF_JZS/BFs is declared
  # inside metrics_arg, it will replace the calibrated one.

  stopifnot(is.logical(calibrated), length(calibrated) == 1)
  stopifnot(is.numeric(alpha), length(alpha) == 1, is.finite(alpha),
            alpha > 0, alpha < 0.5)
  alpha_level <- alpha    # one-sided level of the Family 1 metrics

  stopifnot(is.data.frame(df), is.character(metric))
  stopifnot(all(c("es_o", "se_o", "es_r", "se_r", "study_id") %in% names(df)))
  stopifnot(metric %in% c("significance_r","significance_agg","sum_p","consistency_cir",
                          "consistency_cio","consistency_pi","small_telescopes",
                          "BF_JZS","BF_meta","BFe","BFr","BFs","pScept","snaphybrid",
                          "equivalence","correspondence","Likratio","minimum_effect"))


  # One-sided p-values in the direction of the original claim (Fisher-z) ----
  # Used by the metrics that combine p-values (significance criterion,
  # calibrated meta-analysis, sum of p-values, sceptical p-value).
  #  * Evidence: the reported two-sided p-values (orig_p_value, repli_p_value;
  #    exact values only) when supplied, otherwise es / se.
  #  * Direction of the replication: same_direction when supplied, otherwise
  #    the signs of the estimates.
  # Returns z-values oriented so that the original claim is positive, and the
  # corresponding one-sided p-values. Without the additional columns, this
  # equals the previous computation from es / se.
  oriented_p <- function(df) {
    z_o <- abs(df$es_o / df$se_o)
    z_r <- abs(df$es_r / df$se_r)
    # reported p-values only for pairs with both estimates, so that all metrics
    # are computed for the same claims
    ok <- is.finite(z_o) & is.finite(z_r)
    if ("orig_p_value" %in% names(df)) {
      i <- ok & !is.na(df$orig_p_value)
      z_o[i] <- stats::qnorm(df$orig_p_value[i] / 2, lower.tail = FALSE)
    }
    if ("repli_p_value" %in% names(df)) {
      i <- ok & !is.na(df$repli_p_value)
      z_r[i] <- stats::qnorm(df$repli_p_value[i] / 2, lower.tail = FALSE)
    }
    p_r_own <- stats::pnorm(z_r, lower.tail = FALSE)   # in the direction of the replication estimate
    same_dir <- sign(df$es_o) == sign(df$es_r)
    if ("same_direction" %in% names(df)) {
      same_dir <- ifelse(is.na(df$same_direction), same_dir, df$same_direction)
    }
    z_r <- ifelse(same_dir, z_r, -z_r)
    list(z_o = z_o, z_r = z_r,
         p_o = stats::pnorm(z_o, lower.tail = FALSE),
         p_r = stats::pnorm(z_r, lower.tail = FALSE),
         p_r_own = p_r_own, same_dir = same_dir)
  }

  # SCORE pattern flag as an additional success requirement, when supplied
  # (missing values do not block success). It can only remove successes, so
  # the overall T1E rate stays at most alpha^2.
  pattern_gate <- function(df) {
    if (!"repli_pattern_criteria_met" %in% names(df)) return(rep(TRUE, nrow(df)))
    ifelse(is.na(df$repli_pattern_criteria_met), TRUE, df$repli_pattern_criteria_met)
  }

  # Significance criterion (vote counting) ----
  assess_rep_significance <- function(df, alpha = alpha_level) {

    # Sidedness follows the analysis scale:
    #  fisherz -> ONE-SIDED p in the direction of the estimate, compared to
    #  alpha (original implementation).
    #  cor     -> TWO-SIDED p compared to 2*alpha (numerically the same
    #  decision when combined with the same-direction requirement below)
    #  Setting scale=cor and including the additional arguments (DoF, effective sample size, etc.)
    #  exactly replicates results at the same level of the original SCORE paper.
    #  fisherz: one-sided p-values from oriented_p() (reported p-values and
    #  same_direction when supplied, otherwise es / se), plus the pattern flag
    #  when supplied
    if (unique(df$scale) == "cor") {
      p_o <- 2 * pnorm(abs(df$es_o / df$se_o), lower.tail = FALSE)
      p_r <- 2 * pnorm(abs(df$es_r / df$se_r), lower.tail = FALSE)
      alpha <- 2 * alpha
      same_dir <- sign(df$es_o) == sign(df$es_r)
      gate <- rep(TRUE, nrow(df))
    } else {
      op <- oriented_p(df)
      p_o <- op$p_o
      p_r <- op$p_r_own                 # one-sided p in the direction of the replication estimate
      same_dir <- op$same_dir
      gate <- pattern_gate(df)
    }

    ### Success criterion

    outcome_detailed <- ifelse(p_o >= alpha,
                               "original effect is not significant",
                               ifelse(p_r >= alpha,
                                      "replication effect is not significant",
                                      ifelse(!same_dir,
                                             "direction of effect is inconsistent",
                                             ifelse(!gate,
                                                    "replication effect is significant, but the pattern is not met",
                                                    "replication effect is significant and same direction of original"))))

    outcome <- ifelse(outcome_detailed == "original effect is not significant", "OS not significant",
                      ifelse(outcome_detailed == "replication effect is significant and same direction of original",
                             "Replication success", "Replication failure"))

    data.frame(outcome=outcome, outcome_detailed=outcome_detailed, quantitative_value = p_r)
  }


  # Fixed Meta-Analysis (aggregating one original and one replication) ----
  assess_agg_significance <- function(df, alpha = alpha_level) {

    if (calibrated) {
      # CALIBRATED: unweighted Stouffer, exact overall T1E at alpha^2 ----

      if (unique(df$scale) == "fisherz") {
        # z-values oriented to the original claim (reported p-values and
        # same_direction when supplied, otherwise es / se)
        op   <- oriented_p(df)
        z_o  <- op$z_o
        z_r  <- op$z_r
        gate <- pattern_gate(df)
      } else {
        z_o <- df$es_o / df$se_o
        z_r <- df$es_r / df$se_r
        flip <- ifelse(!is.na(z_o) & z_o < 0, -1, 1)   # orient to the original's direction
        z_o <- z_o * flip
        z_r <- z_r * flip
        gate <- rep(TRUE, nrow(df))
      }

      z_sum  <- z_o + z_r
      thresh <- calibrate_threshold("significance_agg", alpha = alpha)
      p_meta <- stats::pnorm(z_sum / sqrt(2), lower.tail = FALSE)   # combined one-sided p

      success <- !is.na(z_sum) & is.finite(z_sum) & (z_sum >= thresh) & gate

      outcome_detailed <- ifelse(is.na(z_sum) | !is.finite(z_sum), NA_character_,
                                 ifelse(success,
                                        "Aggregated effect exceeds the success threshold",
                                        "Aggregated effect does not exceed the success threshold"))
      outcome <- ifelse(is.na(z_sum) | !is.finite(z_sum), NA_character_,
                        ifelse(success, "Replication success", "Replication failure"))

      return(data.frame(outcome = outcome, outcome_detailed = outcome_detailed,
                        quantitative_value = p_meta))
    }


    # Calculates the pooled effect estimate and the (1 - 2 alpha) CI
    # (95% CI for the default alpha = 0.025)
    effect_comparison <- compare_effectsizes(df$es_o, df$se_o, df$es_r, df$se_r,
                                             level = (1 - 2*alpha)*100)

    ### Calculate p-value for this metric
    p_meta <- effect_comparison$pval

    # Success criterion:
    #   sig = sign(lb_b) == sign(ub_b)      (pooled CI excludes zero)
    #   dir = sign(orig_conv_r) == sign(lb_b)  (direction matches the original)
    # NOTE: this needs signed es_o/ci_ub if you're feeding in abs()
    # correlations this collapses to the old positive-only check, which is
    # wrong for originally-negative effects.
    sig <- sign(effect_comparison$ci_lb) == sign(effect_comparison$ci_ub)
    dir <- sign(df$es_o) == sign(effect_comparison$ci_lb)

    outcome_detailed <- ifelse(is.na(sig) | is.na(dir), NA_character_,
                               ifelse(sig & dir, "Aggregated effect is larger than 0 and in the same direction",
                                      ifelse(sig & !dir, "Aggregated effect is larger than 0 but in the opposite direction",
                                             ifelse(!sig, "Aggregated effect is not larger than 0", NA_character_))))

    outcome <- ifelse(outcome_detailed == "Aggregated effect is larger than 0 and in the same direction",
                      "Replication success",
                      ifelse(outcome_detailed %in% c("Aggregated effect is larger than 0 but in the opposite direction",
                                                     "Aggregated effect is not larger than 0"),
                             "Replication failure", NA_character_))

    data.frame(outcome=outcome, outcome_detailed=outcome_detailed, quantitative_value = p_meta)

  }

  # Sum of $p$-values ----

  assess_sum_p <- function(df, alpha = alpha_level,
                           threshold = calibrate_threshold("sum_p", alpha = alpha)){

    if (unique(df$scale) == "fisherz") {
      # One-sided p-values in the direction of the original claim (reported
      # p-values and same_direction when supplied, otherwise es / se). A
      # replication in the opposite direction has p_r > 0.5 and cannot succeed;
      # the pattern flag is an additional requirement when supplied.
      op <- oriented_p(df)
      sum_p <- op$p_o + op$p_r
      pattern_met <- pattern_gate(df)
    } else {
      if ("orig_p_value"  %in% names(df)){
        p_o <- df$orig_p_value
      } else {
        z_o <- df$es_o / df$se_o
        p_o <- rep(NA_real_, length(z_o))
        p_o[is.finite(z_o)] <- ReplicationSuccess::z2p(z_o[is.finite(z_o)])
      }

      if ("repli_p_value"  %in% names(df)){
        p_r <- df$repli_p_value
      } else {
        z_r <- df$es_r / df$se_r
        p_r <- rep(NA_real_, length(z_r))
        p_r[is.finite(z_r)] <- ReplicationSuccess::z2p(z_r[is.finite(z_r)])
      }

      sum_p <- (p_o / 2) + (p_r / 2)

      ### the sum-of-p test is direction-agnostic, so AND with the qualitative
      ### `sum_p <- both_p & repli_pattern_criteria_met`
      pattern_met <- if ("repli_pattern_criteria_met" %in% names(df)) {
        df$repli_pattern_criteria_met
      } else {
        sign(df$es_o) == sign(df$es_r)
      }
    }

    both_p <- sum_p < threshold

    success <- !is.na(both_p) & !is.na(pattern_met) & both_p & pattern_met

    outcome_detailed <- ifelse(is.na(both_p) | is.na(pattern_met), NA_character_,
                               ifelse(success,
                                      "Sum of one-tailed p-values is below threshold and pattern matches",
                                      ifelse(both_p & !pattern_met,
                                             "Sum of one-tailed p-values is below threshold but pattern does not match",
                                             "Sum of one-tailed p-values is not below threshold")))

    outcome <- ifelse(is.na(both_p) | is.na(pattern_met), NA_character_,
                      ifelse(success, "Replication success", "Replication failure"))

    data.frame(outcome = outcome, outcome_detailed = outcome_detailed,
               quantitative_value = sum_p)
  }

  # Original effect in Replication Confidence interval metric ----
  # `alpha`: one-sided level, (1 - 2 alpha) interval. Reported 95% limits are
  # used as they are at alpha = 0.025 and rescaled otherwise (rescale_ci()).
  assess_or_consistency_with_rep <- function(df, alpha = 0.025) {

    if (unique(df$scale) == "cor" && all(c("repli_conv_r_lb","repli_conv_r_ub") %in% names(df))) {

      ci <- rescale_ci(df$repli_conv_r_lb, df$repli_conv_r_ub, df$es_r, alpha = alpha)
      ci <- fill_ci(ci, df$es_r, df$se_r, alpha)   # rows without reported limits
      outcome <- ifelse(df$es_o > ci$lower & df$es_o < ci$upper,
                        "Replication success", "Replication failure")

      data.frame(outcome=outcome, outcome_detailed=outcome, quantitative_value=NA_real_)

    } else if (unique(df$scale) == "fisherz" && all(c("es_lb_r","es_ub_r") %in% names(df))) {
      # Containment against supplied bounds on the z scale (atanh of the study's
      # reported CI). atanh is monotonic, so decision-identical to raw-r.
      ci <- rescale_ci(df$es_lb_r, df$es_ub_r, df$es_r, alpha = alpha, on_z = TRUE)
      ci <- fill_ci(ci, df$es_r, df$se_r, alpha)   # rows without reported limits
      outcome <- ifelse(is.na(df$es_o) | is.na(ci$lower) | is.na(ci$upper), NA_character_,
                        ifelse(df$es_o > ci$lower & df$es_o < ci$upper,
                               "Replication success", "Replication failure"))
      data.frame(outcome=outcome, outcome_detailed=outcome, quantitative_value=NA_real_)

    } else {

      p_cir <- pnorm((df$es_r - df$es_o) / df$se_r, mean = 0, sd = 1, lower.tail = FALSE)
      cis <- compute_ci(df$es_r, df$se_r, alpha = alpha)

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
  }

  # Replication in Original Confidence interval metric ----
  assess_rep_consistency_with_or <- function(df, alpha = 0.025) {

    if (unique(df$scale) == "cor" && all(c("orig_conv_r_lb","orig_conv_r_ub") %in% names(df))) {
      # --- precomputed bounds (rescaled when alpha != 0.025)
      ci <- rescale_ci(df$orig_conv_r_lb, df$orig_conv_r_ub, df$es_o, alpha = alpha)
      ci <- fill_ci(ci, df$es_o, df$se_o, alpha)   # rows without reported limits
      outcome <- ifelse(df$es_r > ci$lower & df$es_r < ci$upper,
                        "Replication success", "Replication failure")
      data.frame(outcome=outcome, outcome_detailed=outcome, quantitative_value=NA_real_)

    } else if (unique(df$scale) == "fisherz" && all(c("es_lb_o","es_ub_o") %in% names(df))) {
      # Same monotonicity argument: decision-identical to raw-r.
      ci <- rescale_ci(df$es_lb_o, df$es_ub_o, df$es_o, alpha = alpha, on_z = TRUE)
      ci <- fill_ci(ci, df$es_o, df$se_o, alpha)   # rows without reported limits
      outcome <- ifelse(is.na(df$es_r) | is.na(ci$lower) | is.na(ci$upper), NA_character_,
                        ifelse(df$es_r > ci$lower & df$es_r < ci$upper,
                               "Replication success", "Replication failure"))
      data.frame(outcome=outcome, outcome_detailed=outcome, quantitative_value=NA_real_)

    } else {
      # --- existing Fisher-z / self-computed-CI path ---
      p_cio <- pnorm((df$es_o - df$es_r) / df$se_o, mean = 0, sd = 1, lower.tail = FALSE)
      cis <- compute_ci(df$es_o, df$se_o, alpha = alpha)
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
  }

  # Replication in Original Prediction Interval ----
  assess_rep_consistency_with_pi <- function(df, alpha = 0.025) {

    ### Calculate p-value for this metric
    p_PI <- pnorm((df$es_o - df$es_r) / sqrt(df$se_o^2 + df$se_r^2), mean = 0, sd = 1, lower.tail = FALSE)

    ### Calculate the (1 - 2 alpha) prediction interval
    pred_intervals <- compute_prediction_interval(df$es_o, df$se_o, df$se_r, alpha = alpha)

    # On the correlation scale, cap the bounds to [-1, 1].
    # Not applicable on the Fisher-z scale,
    # where the effect size isn't bounded.
    if (unique(df$scale) == "cor") {
      pred_intervals$lower <- pmax(pred_intervals$lower, -1)
      pred_intervals$upper <- pmin(pred_intervals$upper,  1)
    }

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

  # `alpha`: level of the test of the replication against the benchmark effect.
  # `alpha_benchmark`: one-sided level used to define the benchmark (the effect
  # the original had `power` to detect), on every path.
  # The test is always carried out in the direction of the original estimate:
  # the replication estimate (and reported limits) are mirrored when the
  # original estimate is negative. Rows with reported limits compare them with
  # a pwr.r.test() benchmark (rescaled when alpha != 0.025); rows without
  # reported limits use the probabilistic z-test path and its benchmark.
  assess_small_telescopes <- function(df, power=1/3, alpha=0.025, alpha_benchmark = alpha) {

    n1_o <- if ("n1_o" %in% names(df)) df$n1_o else rep(NA_real_, nrow(df))
    n2_o <- if ("n2_o" %in% names(df)) df$n2_o else rep(NA_real_, nrow(df))
    scale_here <- unique(df$scale)

    # +1 when the original estimate is positive (or zero), -1 when negative,
    # NA when it is missing (the pair cannot be assessed)
    flip <- ifelse(is.na(df$es_o), NA_real_, ifelse(df$es_o < 0, -1, 1))

    # Probabilistic path: benchmark from a one-sided z power calculation on
    # the Fisher-z scale (back-transformed to r when scale = "cor")
    bench_prob <- mapply(function(effectsize, n1, n2, n) {
      if (effectsize %in% c("SMD", "LOGOR")) {
        if (is.na(n1) || is.na(n2)) return(NA_real_)
        d33 <- find_z_for_power(n1 = n1, n2 = n2, power = power, alpha = alpha_benchmark)
        r_to_zr(d_to_r(d33, n1 = n1, n2 = n2)[["r"]])[[1]]
      } else {
        if (is.na(n) || n < 4) return(NA_real_)
        power.z.test(power = power, n = n-3, type = "one.sample",
                     sig.level = alpha_benchmark, alternative = "one.sided")$delta
      }
    }, effectsize=df$effectsize, n1=n1_o, n2=n2_o, n=df$n_o)
    bench_prob <- as.numeric(bench_prob)
    if (scale_here == "cor") bench_prob <- tanh(bench_prob)

    p_telescopes <- pnorm((flip * df$es_r - bench_prob) / df$se_r, mean=0, sd=1, lower.tail=TRUE)
    hit <- p_telescopes >= alpha
    quant_val <- p_telescopes

    # Rows with reported limits (on either side)
    reported_cor <- scale_here == "cor" && all(df$effectsize == "COR") &&
      all(c("repli_conv_r_lb", "repli_conv_r_ub") %in% names(df))
    reported_z <- scale_here == "fisherz" &&
      all(c("es_lb_r", "es_ub_r") %in% names(df)) &&
      all(df$effectsize %in% c("COR", "ZCOR"))

    if (reported_cor || reported_z) {
      bench_r <- function(n) {
        if (is.na(n) || !is.finite(n) || n < 4) return(NA_real_)
        tryCatch(pwr::pwr.r.test(n = n, power = power, sig.level = alpha_benchmark,
                                 alternative = "greater")$r,
                 error = function(e) NA_real_)
      }
      lb <- if (reported_cor) df$repli_conv_r_lb else df$es_lb_r
      ub <- if (reported_cor) df$repli_conv_r_ub else df$es_ub_r
      has_rep <- !is.na(lb) | !is.na(ub)
      if (any(has_rep)) {
        small_effect <- vapply(df$n_o, bench_r, numeric(1))
        if (reported_z) small_effect <- atanh(small_effect)
        ci <- rescale_ci(lb, ub, df$es_r, alpha = alpha, on_z = reported_z)
        ci <- fill_ci(ci, df$es_r, df$se_r, alpha)   # a single missing side
        # upper limit of the replication CI in the direction of the original
        upper_oriented <- ifelse(flip == 1, ci$upper, -ci$lower)
        hit[has_rep] <- (upper_oriented > small_effect)[has_rep]
        quant_val[has_rep] <- small_effect[has_rep]
      }
    }

    outcome_detailed <- ifelse(is.na(hit), NA_character_,
                               ifelse(hit,
                                      "original study had >= 33% power to detect replication effect",
                                      "original study had < 33% power to detect replication effect"))
    outcome <- ifelse(is.na(hit), NA_character_,
                      ifelse(hit, "Replication success", "Replication failure"))

    data.frame(outcome=outcome, outcome_detailed=outcome_detailed, quantitative_value = quant_val)
  }

  # JZS Bayes Factor ----
  assess_BFJZS <- function(df, gamma = 1){

    # Calibrated threshold with overall T1E rate at most alpha^2 (for a
    # pre-specified direction) for every standard error (1.60 for
    # alpha = 0.025), unless gamma is given explicitly (metric_args)
    if (calibrated && missing(gamma)) gamma <- calibrate_threshold("BF_JZS", alpha = alpha_level)

    # log JZS Bayes factor (Cauchy(0, 1) prior on the effect); the integral is
    # computed by .bfjzs_log() over es +/- 40 se, which stays accurate for
    # very small standard errors
    compute_bfjzs <- function(es, se) {

      # Handle missing or invalid inputs
      if (is.na(es) || is.na(se) ||
          !is.finite(es) || !is.finite(se) ||
          se <= 0) {
        return(NA_real_)
      }

      tryCatch(.bfjzs_log(es, se), error = function(e) NA_real_)
    }

    BFJZS_ori <- mapply(compute_bfjzs, es = df$es_o, se = df$se_o)
    BFJZS_rep <- mapply(compute_bfjzs, es = df$es_r, se = df$se_r)

    # Direction check (one-sided): success requires the replication estimate to
    # have the direction of the original; evidence for an effect in the
    # opposite direction is a failure
    same_dir <- sign(df$es_r) == sign(df$es_o)

    ### Success criterion

    lg <- log(gamma)
    ok <- !is.na(BFJZS_ori) & !is.na(BFJZS_rep) & !is.na(same_dir)

    outcome_detailed <- ifelse(!ok, NA_character_,
                               ifelse(BFJZS_ori < -lg, "Evidence for no original effect",
                                      ifelse(BFJZS_ori <= lg, "Original null effect",
                                             ifelse(BFJZS_rep > lg & same_dir, "Evidence for replication success",
                                                    ifelse(BFJZS_rep > lg, "Evidence for an effect opposite to the original",
                                                           ifelse(BFJZS_rep < -lg, "Evidence for replication failure",
                                                                  "Inconclusive replication"))))))

    outcome <- ifelse(!ok, NA_character_,
                      ifelse(outcome_detailed == "Evidence for replication success", "Replication success",
                             ifelse(outcome_detailed %in% c("Evidence for an effect opposite to the original",
                                                            "Evidence for replication failure"),
                                    "Replication failure", "Inconclusive")))
    data.frame(outcome=outcome, outcome_detailed=outcome_detailed, quantitative_value = exp(BFJZS_rep))
  }

  # Fixed-effect meta-analytic Bayes factor ----
  # Rouder and Morey (2011); Verhagen and Wagenmakers (2014). H0: common effect
  # theta = 0 versus H1: theta ~ Cauchy(0, 1), the same effect in the original
  # and replication studies. With normal likelihoods, the product of the two
  # likelihoods factorises into a term that does not depend on theta and the
  # likelihood of the inverse-variance pooled estimate, so the Bayes factor is
  # the JZS Bayes factor of the pooled estimate (computed by .bfjzs_log(), see helper.R).
  assess_BFmeta <- function(df, gamma = 1){

    # Calibrated threshold with overall T1E rate at most alpha^2 for every
    # standard error of the pooled estimate, unless gamma is given explicitly (metric_args)
    if (calibrated && missing(gamma)) gamma <- calibrate_threshold("BF_meta", alpha = alpha_level)

    ok   <- is.finite(df$es_o) & is.finite(df$se_o) & df$se_o > 0 &
      is.finite(df$es_r) & is.finite(df$se_r) & df$se_r > 0
    w_o  <- 1 / df$se_o^2
    w_r  <- 1 / df$se_r^2
    es_p <- (w_o * df$es_o + w_r * df$es_r) / (w_o + w_r)   # pooled estimate
    se_p <- sqrt(1 / (w_o + w_r))                           # its standard error

    logBF <- rep(NA_real_, nrow(df))
    if (any(ok)) {
      logBF[ok] <- mapply(function(es, se) tryCatch(.bfjzs_log(es, se), error = function(e) NA_real_),
                          es_p[ok], se_p[ok], USE.NAMES = FALSE)
    }

    # One-sided: success requires the pooled estimate to have the direction of
    # the original; evidence for a pooled effect in the opposite direction is a
    # failure
    same_dir <- sign(es_p) == sign(df$es_o)

    ### Success criterion

    outcome_detailed <- ifelse(is.na(logBF), NA_character_,
                               ifelse(logBF > log(gamma) & same_dir,
                                      "Evidence for a pooled effect in the direction of the original",
                                      ifelse(logBF > log(gamma),
                                             "Evidence for a pooled effect opposite to the original",
                                             ifelse(logBF < log(1/gamma),
                                                    "Evidence for no pooled effect",
                                                    "Inconclusive evidence"))))

    outcome <- ifelse(is.na(logBF), NA_character_,
                      ifelse(logBF > log(gamma) & same_dir, "Replication success",
                             ifelse(logBF > log(gamma) | logBF < log(1/gamma), "Replication failure",
                                    "Inconclusive")))

    data.frame(outcome = outcome, outcome_detailed = outcome_detailed,
               quantitative_value = exp(logBF))
  }

  # Equality of effect Bayes Factor ----
  assess_BFe <- function(df, gamma = 1, tau=0.2){

    compute_BFe <- function(to, so, tr, sr, tau){

      if(any(is.na(c(to, so, tr, sr))) ||
         any(!is.finite(c(to, so, tr, sr))) ||
         so <= 0 || sr <= 0){
        return(NA_real_)
      }

      tryCatch(
        BayesRep::BFe(to = to, so = so, tr = tr, sr = sr, tau = tau),
        error = function(e) NA_real_
      )
    }

    BFe <- mapply(compute_BFe, to = df$es_o, so = df$se_o,
                  tr = df$es_r, sr = df$se_r, tau = tau)

    # ## compute sceptical BF and invert it, so that values >1 indicate evidence consistent with original study (replication success)
    # BFe <- BayesRep::BFe(to=df$es_o, so=df$se_o, tr=df$es_r, sr=df$se_r, tau=tau)

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

  assess_BFr <- function(df, gamma = 1){

    compute_BFr <- function(to, so, tr, sr){
      if(any(is.na(c(to, so, tr, sr))) ||
         any(!is.finite(c(to, so, tr, sr))) ||
         so <= 0 || sr <= 0){
        return(NA_real_)
      }

      tryCatch(
        BayesRep::BFr(to = to, so = so, tr = tr, sr = sr),
        error = function(e) NA_real_
      )
    }

    ## compute sceptical BF and invert it, so that values >1 indicate evidence consistent with original study (replication success)

    BFr <- 1 / mapply(compute_BFr, to = df$es_o, so = df$se_o,
                      tr = df$es_r, sr = df$se_r)

    ### Success criterion

    # `bayes_rep <- repli_pattern_criteria_met & r_bf < 1`
    pattern_met <- if ("repli_pattern_criteria_met" %in% names(df)) {
      df$repli_pattern_criteria_met
    } else {
      sign(df$es_o) == sign(df$es_r)
    }

    success <- !is.na(BFr) & !is.na(pattern_met) & (BFr > gamma) & pattern_met

    outcome_detailed <- ifelse(is.na(BFr) | is.na(pattern_met), NA_character_,
                               ifelse(success,       "Evidence for replication success",
                                      ifelse(BFr >= 1,           "Weak evidence for replication success",
                                             ifelse(BFr >= 1/gamma,     "Weak evidence for replication failure",
                                                    ifelse(BFr <= 1/gamma,     "Evidence for replication failure",
                                                           NA_character_)))))

    outcome <- ifelse(is.na(BFr) | is.na(pattern_met), NA_character_,
                      ifelse(success,    "Replication success",
                             ifelse(BFr <= 1/gamma, "Replication failure",
                                    "Inconclusive")))

    data.frame(outcome=outcome, outcome_detailed=outcome_detailed, quantitative_value = BFr)
  }

  # Sceptical Bayes Factor ----

  assess_BFs <- function(df, gamma = 1){

    # Calibrated threshold: the T1E rate of the sceptical Bayes factor depends
    # only on the variance ratio c = se_o^2 / se_r^2, so the threshold is
    # computed for every pair from its own c, which gives an overall T1E rate
    # of exactly alpha^2 for every pair (2.36 for alpha = 0.025 and c = 1).
    # An explicitly given gamma (metric_args) replaces the calibrated one.
    # Pairs without a sceptical Bayes factor (NaN) count as non-successes.
    if (calibrated && missing(gamma)) {
      c_vr  <- df$se_o^2 / df$se_r^2
      ok    <- is.finite(c_vr) & c_vr > 0
      gamma <- rep(NA_real_, nrow(df))
      if (any(ok)) gamma[ok] <- calibrate_threshold("BFs", alpha = alpha_level, c = c_vr[ok])
    }

    compute_BFs <- function(to, so, tr, sr){
      if(any(is.na(c(to, so, tr, sr))) ||
         any(!is.finite(c(to, so, tr, sr))) ||
         so <= 0 || sr <= 0){
        return(NA_real_)
      }

      tryCatch(
        BayesRep::BFs(to = to, so = so, tr = tr, sr = sr),
        error = function(e) NA_real_
      )
    }

    ## compute sceptical BF and invert it, so that values >1 indicate evidence consistent with original study (replication success)
    BFs <- 1 / mapply(compute_BFs, to=df$es_o, so=df$se_o,
                      tr=df$es_r, sr=df$se_r)

    ### Success criterion

    same_dir <- sign(df$es_r) == sign(df$es_o)

    outcome_detailed <- ifelse(is.nan(BFs),
                               "Results of replication is too inconclusive to establish replication success at any level",
                               ifelse(BFs > gamma & same_dir, "Evidence for replication success",
                                      ifelse(BFs > gamma,            "Replication estimate in the opposite direction",
                                             ifelse(BFs >= 1,               "Weak evidence for replication success",
                                                    ifelse(BFs >= 1/gamma,         "Weak evidence for replication failure",
                                                           ifelse(BFs < 1/gamma,          "Evidence for replication failure",
                                                                  NA_character_))))))

    outcome <- ifelse(is.nan(BFs), "Inconclusive",
                      ifelse(BFs > gamma & same_dir, "Replication success",
                             ifelse(BFs > gamma | BFs <= 1/gamma, "Replication failure",
                                    "Inconclusive")))

    data.frame(outcome=outcome, outcome_detailed=outcome_detailed, quantitative_value = BFs)
  }

  # Sceptical p-value ----

  assess_pScept <- function(df, alpha = alpha_level){

    # Sidedness follows the analysis scale:
    #  scale == "cor"     -> TWO-SIDED: z's from the reported p-values via
    #    p2z(); success gated by the qualitative pattern flag (two.sided
    #    pSceptical doesn't check direction itself).
    # scale == "fisherz" -> ONE-SIDED: z's from the effect size / SE ratio;
    #    no pattern (one.sided pSceptical handles direction already).
    #
    # pSceptical(type = "golden" / "controlled") returns a RECALIBRATED sceptical
    # p-value, which is compared directly with the significance level: alpha
    # (one-sided, fisher z) or 2 * alpha (two-sided, cor scale).
    # levelSceptical() gives the threshold for the NOMINAL p-value and must not
    # be applied to the recalibrated one (that would recalibrate twice).
    score_mode <- unique(df$scale) == "cor"

    # Calibrated: use the "controlled" recalibration, which achieves exact
    # overall T1E control at alpha^2 for every variance ratio c = n_r/n_o;
    # otherwise use the "golden" recalibration.
    ps_type <- if (calibrated) "controlled" else "golden"

    compute_pScept_cor <- function(es_r, se_r, es_o, se_o, n_o, n_r,
                                   p_o = NA_real_, p_r = NA_real_){
      if (any(is.na(c(n_o, n_r))) || n_o <= 0 || n_r <= 0) return(NA_real_)
      p_ok <- function(p) !is.na(p) && is.finite(p) && p > 0 && p <= 1
      if (p_ok(p_o) && p_ok(p_r)) {
        zo <- ReplicationSuccess::p2z(p_o)
        zr <- ReplicationSuccess::p2z(p_r)
      } else {
        return(NA_real_)   # no valid exact p -> unusable -> NA
      }
      tryCatch(
        ReplicationSuccess::pSceptical(zo = zo, zr = zr, c = n_r/n_o,
                                       alternative = "two.sided", type = ps_type),
        error = function(e) NA_real_
      )
    }

    compute_pScept_z <- function(es_r, se_r, es_o, se_o, n_o, n_r){
      if (any(is.na(c(es_r, se_r, es_o, se_o, n_o, n_r))) ||
          any(!is.finite(c(es_r, se_r, es_o, se_o, n_o, n_r))) ||
          se_o <= 0 || se_r <= 0 || n_o <= 0 || n_r <= 0){
        return(NA_real_)
      }
      tryCatch(
        ReplicationSuccess::pSceptical(zo = es_o/se_o, zr = es_r/se_r, c = n_r/n_o,
                                       alternative = "one.sided", type = ps_type),
        error = function(e) NA_real_
      )
    }

    # guard against dropped columns (COR-only data has no n1_o/n2_o/n1_r/n2_r)
    n1_o <- if ("n1_o" %in% names(df)) df$n1_o else rep(NA_real_, nrow(df))
    n2_o <- if ("n2_o" %in% names(df)) df$n2_o else rep(NA_real_, nrow(df))
    n1_r <- if ("n1_r" %in% names(df)) df$n1_r else rep(NA_real_, nrow(df))
    n2_r <- if ("n2_r" %in% names(df)) df$n2_r else rep(NA_real_, nrow(df))

    # resolve n_o/n_r ONCE, vectorized
    n_o_resolved <- ifelse(!is.na(n1_o) & !is.na(n2_o), n1_o + n2_o, df$n_o)
    n_r_resolved <- ifelse(!is.na(n1_r) & !is.na(n2_r), n1_r + n2_r, df$n_r)

    if (score_mode) {
      p_o <- if ("orig_p_value"  %in% names(df)) df$orig_p_value  else rep(NA_real_, nrow(df))
      p_r <- if ("repli_p_value" %in% names(df)) df$repli_p_value else rep(NA_real_, nrow(df))

      pScept <- mapply(compute_pScept_cor,
                       es_o = df$es_o, se_o = df$se_o, es_r = df$es_r, se_r = df$se_r,
                       n_o = n_o_resolved, n_r = n_r_resolved,
                       p_o = p_o, p_r = p_r)

      # two-sided recalibrated p-value -> compare with 2 * alpha
      skep_level <- rep(2 * alpha, nrow(df))

      pattern_met <- if ("repli_pattern_criteria_met" %in% names(df)) {
        df$repli_pattern_criteria_met
      } else {
        sign(df$es_o) == sign(df$es_r)
      }

      success <- !is.na(pScept) & !is.na(pattern_met) & !is.na(skep_level) & (pScept < skep_level) & pattern_met
      na_mask <- is.na(pScept) | is.na(pattern_met)

    } else {
      # z-values oriented to the original claim (reported p-values and
      # same_direction when supplied, otherwise es / se); variance ratio
      # c = n_r / n_o for n-based SEs, otherwise se_o^2 / se_r^2 (SEs from the
      # degrees of freedom or the effective sample size, or supplied SEs)
      op   <- oriented_p(df)
      c_vr <- n_r_resolved / n_o_resolved
      if ("se_src_o" %in% names(df)) {
        # SEs not derived from n (df, effective n, or supplied with ZCOR/SMD/LOGOR
        # input): use the variance ratio of the SEs
        meta <- df$se_src_o %in% c("df", "effective_n", "supplied") |
          df$se_src_r %in% c("df", "effective_n", "supplied")
        c_vr[meta] <- (df$se_o^2 / df$se_r^2)[meta]
      }
      pScept <- mapply(function(zo, zr, cc) {
        if (any(!is.finite(c(zo, zr, cc))) || cc <= 0) return(NA_real_)
        tryCatch(ReplicationSuccess::pSceptical(zo = zo, zr = zr, c = cc,
                                                alternative = "one.sided", type = ps_type),
                 error = function(e) NA_real_)
      }, op$z_o, op$z_r, c_vr)

      # one-sided recalibrated p-value -> compare with alpha
      skep_level <- rep(alpha, nrow(df))

      success <- !is.na(pScept) & !is.na(skep_level) & (pScept <= skep_level) & pattern_gate(df)
      na_mask <- is.na(pScept)
    }

    outcome_detailed <- ifelse(na_mask, NA_character_,
                               ifelse(success,
                                      "Evidence for replication success",
                                      "No evidence for replication success"))
    outcome <- ifelse(na_mask, NA_character_,
                      ifelse(success, "Replication success", "Replication failure"))

    data.frame(outcome=outcome, outcome_detailed=outcome_detailed, quantitative_value = pScept)
  }

  # Snapshot hybrid ----

  assess_snaphybrid <- function(df, threshold=0.75, pr_0 = 0.25, pr_s = 0.25, pr_m = 0.25){
    # Calculate prior for large effect based on other priors
    pr_l <- 1 - (pr_0 + pr_s + pr_m)
    if (is.finite(pr_l) && pr_l < 0 && pr_l > -1e-10) pr_l <- 0   # rounding
    priors <- c(pr_0, pr_s, pr_m, pr_l)
    if (any(!is.finite(priors)) || any(priors < 0) || any(priors > 1))
      stop("snaphybrid: pr_0, pr_s and pr_m must be probabilities with pr_0 + pr_s + pr_m <= 1")

    effect_labels <- c("null effect", "small effect", "medium effect", "large effect")

    base_eligible <- if ("snaphybrid_base_eligible" %in% names(df)) {
      df$snaphybrid_base_eligible
    } else {
      rep(TRUE, nrow(df))
    }
    os_nonsig <- rep(FALSE, nrow(df))
    pSnap   <- rep(NA_real_, nrow(df))
    snap_idx <- rep(NA_integer_, nrow(df))

    for (i in seq_len(nrow(df))) {

      ## Following van Aert (author of snapshot hybrid) The snapshot is applied
      # only when the original is significant (one-sided p < 0.025).
      # That leaves 211 claims
      if (!isTRUE(base_eligible[i])) {
        next
      }

      n_o <- df$n_o[i]
      n_r <- df$n_r[i]

      if (anyNA(c(n_o, n_r)) ||
          any(!is.finite(c(n_o, n_r))) ||
          n_o <= 3 || n_r <= 3) {
        next
      }

      ## Recover raw correlations
      if (df$scale[i] == "fisherz") {
        r_o <- tanh(df$es_o[i])
        r_r <- tanh(df$es_r[i])
      } else if (df$scale[i] == "cor") {
        r_o <- df$es_o[i]
        r_r <- df$es_r[i]
      } else {
        next
      }

      if (anyNA(c(r_o, r_r)) ||
          any(!is.finite(c(r_o, r_r))) ||
          any(abs(c(r_o, r_r)) >= 1)) {
        next
      }

      ## Snapshot uses Fisher-z correlations and raw sample-size variances,
      ## unless the SE comes from degrees of freedom or effective sample size
      es <- atanh(c(r_o, r_r))
      se <- 1 / sqrt(c(n_o, n_r) - 3)
      if (df$scale[i] == "fisherz" && "se_src_o" %in% names(df)) {
        if (!is.na(df$se_src_o[i]) && df$se_src_o[i] %in% c("df", "effective_n")) se[1] <- df$se_o[i]
        if (!is.na(df$se_src_r[i]) && df$se_src_r[i] %in% c("df", "effective_n")) se[2] <- df$se_r[i]
      }

      ## Orient the original positively and reverse the replication with it
      if (es[1] < 0) {
        es <- -es
      }

      ## one-sided p < .025 after orientation (snapshot only when the original is significant)
      p_original <- stats::pnorm(
        es[1] / se[1],
        lower.tail = FALSE
      )

      if (!is.finite(p_original) || p_original >= 0.025) {
        os_nonsig[i] <- TRUE
        next
      }

      vals <- tryCatch(
        hybrid_snap(
          es = es,
          se = se,
          alpha = 0.025,
          priors = priors
        ),
        error = function(e) NULL
      )

      if (is.null(vals) ||
          !is.data.frame(vals) ||
          nrow(vals) != 1L ||
          ncol(vals) != 4L) {
        next
      }

      probs <- as.numeric(vals[1, ])

      if (length(probs) != 4L ||
          anyNA(probs) ||
          any(!is.finite(probs)) ||
          any(probs < 0 | probs > 1) ||
          abs(sum(probs) - 1) > 1e-6) {
        next
      }

      idx <- which.max(probs)

      pSnap[i]    <- probs[idx]
      snap_idx[i] <- idx
    }

    evaluated       <- !is.na(pSnap) & !is.na(snap_idx)
    above_threshold <- evaluated & pSnap > threshold
    nonnull_selected <- above_threshold & snap_idx %in% 2:4

    outcome_detailed <- rep(NA_character_, nrow(df))

    outcome_detailed[evaluated & !above_threshold] <-
      "Inconclusive evidence"

    outcome_detailed[above_threshold] <- paste(
      "Evidence for",
      effect_labels[snap_idx[above_threshold]]
    )

    ## Collapse null and inconclusive results into replication failure
    outcome <- rep(NA_character_, nrow(df))
    outcome[evaluated] <- "Replication failure"
    outcome[nonnull_selected] <- "Replication success"
    outcome[os_nonsig]          <- "OS not significant"
    outcome_detailed[os_nonsig] <- "original effect is not significant"

    data.frame(
      outcome = outcome,
      outcome_detailed = outcome_detailed,
      quantitative_value = pSnap
    )
  }

  # snap_results <- data.frame(
  #   pSnap = rep(NA_real_, nrow(df)),
  #   snap_idx = rep(NA_integer_, nrow(df))
  # )
  #
  # for(i in seq_len(nrow(df))){
  #
  #   es_o <- df$es_o[i]
  #   es_r <- df$es_r[i]
  #   se_o <- df$se_o[i]
  #   se_r <- df$se_r[i]
  #
  #   ## Skip missing values
  #   if(any(is.na(c(es_o, es_r, se_o, se_r))) ||
  #      any(!is.finite(c(es_o, es_r, se_o, se_r))) ||
  #      se_o <= 0 || se_r <= 0){
  #     next
  #   }
  #
  #   vals <- tryCatch(
  #     hybrid_snap(
  #       es = c(es_o, es_r),
  #       se = c(se_o, se_r)
  #     ),
  #     error = function(e) NULL
  #   )
  #
  #   ## hybrid_snap failed
  #   if(is.null(vals))
  #     next
  #
  #   ## Check output format
  #   if(!is.data.frame(vals) ||
  #      nrow(vals) != 1 ||
  #      ncol(vals) != 4)
  #     next
  #
  #   probs <- as.numeric(vals[1, ])
  #
  #   ## Invalid probabilities
  #   if(anyNA(probs) ||
  #      length(probs) != 4)
  #     next
  #
  #   idx <- which.max(probs)
  #
  #   snap_results$pSnap[i] <- probs[idx]
  #   snap_results$snap_idx[i] <- idx
  # }
  #
  # cat("Rows in df:", nrow(df), "\n")
  # cat("Rows in snap_results:", nrow(snap_results), "\n")
  #
  # post_Snap <- snap_results$pSnap
  # snap_idx  <- snap_results$snap_idx
  #
  # ## Rows for which snapshot hybrid was successfully calculated
  # evaluated <- !is.na(post_Snap) & !is.na(snap_idx)
  #
  #
  # above_threshold <- evaluated & post_Snap > threshold
  #
  # ## Selected model
  # null_selected <- above_threshold & snap_idx == 1L
  # nonnull_selected <- above_threshold & snap_idx %in% 2:4
  #
  # ### Detailed classification
  #
  # outcome_detailed <- rep(NA_character_, nrow(df))
  #
  # ## Evaluated rows without any PMP > threshold
  # outcome_detailed[evaluated & !above_threshold] <-
  #   "Inconclusive evidence"
  #
  # ## At threshold = 0.75, at most one model can satisfy this condition
  # outcome_detailed[above_threshold] <- paste(
  #   "Evidence for",
  #   effect_labels[snap_idx[above_threshold]]
  # )
  #
  # ### Binary replication outcome
  #
  # outcome <- rep(NA_character_, nrow(df))
  #
  # ## Null evidence and inconclusive evidence are both failures
  # outcome[evaluated] <- "Replication failure"
  #
  # ## Small, medium, or large effect evidence is success
  # outcome[nonnull_selected] <- "Replication success"
  #
  # data.frame(outcome = outcome,
  #   outcome_detailed = outcome_detailed,
  #   quantitative_value = post_Snap)}

  # Likelihood ratio approach ----

  # assess_Likratio <- function(df, lr = 8){
  #   Dtie_res <- Dtie(n_o = n_o, lr = 8, Correction = Qc)
  #
  #   #
  #   LRae_0 <- LR.ae_0(es_r = es_r, Dtie = Dtie_res, n_r = n_r)
  #
  #   ### Success criterion
  #
  #   outcome_detailed <- ifelse(LRae_0 >= lr,
  #                              "Evidence for replication success",
  #                              ifelse(LRae_0 > 1/lr && LRae_0 < lr,
  #                                     "Inconclusive evidence",
  #                                     ifelse(LRae_0 <= 1/lr,
  #                                            "Evidence for replication failure",
  #                                            NA_character_)))
  #
  #
  #   outcome <- ifelse(LRae_0 >= lr,
  #                     "Evidence for replication success",
  #                     ifelse(LRae_0 > 1/lr && LRae_0 < lr,
  #                            "Inconclusive evidence",
  #                            ifelse(LRae_0 <= 1/lr,
  #                                   "Evidence for replication failure",
  #                                   NA_character_)))
  #
  #   data.frame(outcome=outcome, outcome_detailed=outcome_detailed, quantitative_value = LRae_0)
  #
  # }

  # `lr`: evidence threshold. `lr_dtie`: likelihood ratio that defines the
  # theoretically interesting effect Dtie (by default the same as lr).
  assess_Likratio <- function(df, lr = 8, lr_dtie = lr){

    scale_here <- unique(df$scale)
    to_r    <- function(es) if (identical(scale_here, "fisherz")) tanh(es) else es
    r_to_d  <- function(r)  2 * r / sqrt(1 - r^2)

    compute_Likratio <- function(es_o, es_r, n_o, n_r, lr_dtie) {
      if (any(is.na(c(es_o, es_r, n_o, n_r))) ||
          any(!is.finite(c(es_o, es_r, n_o, n_r))) ||
          n_o <= 2 || n_r <= 2) return(NA_real_)

      r_o <- to_r(es_o)
      r_r <- to_r(es_r)
      if (any(is.na(c(r_o, r_r))) || abs(r_o) >= 1 || abs(r_r) >= 1) return(NA_real_)

      # mirror so the replication is expressed in the original's direction
      if (r_o < 0) r_r <- -r_r

      d_r <- r_to_d(r_r)
      if (!is.finite(d_r)) return(NA_real_)

      # Dtie() ends in sqrt() and returns NA for very small n_o
      # (the variance ratio drops below 1). is.na() handles this potential error
      Dtie_res <- suppressWarnings(
        tryCatch(Dtie(n_o = n_o, lr = lr_dtie), error = function(e) NA_real_)
      )
      if (is.na(Dtie_res) || !is.finite(Dtie_res)) return(NA_real_)

      out <- suppressWarnings(
        tryCatch(LR.ae_0(es_r = d_r, Dtie = Dtie_res, n_r = n_r),
                 error = function(e) NA_real_)
      )
      if (!is.finite(out)) return(NA_real_)
      out
    }

    LRae_0 <- mapply(compute_Likratio,
                     es_o = df$es_o, es_r = df$es_r,
                     n_o = df$n_o, n_r = df$n_r,
                     MoreArgs = list(lr_dtie = lr_dtie))

    ### Success criterion --
    outcome_detailed <- ifelse(is.na(LRae_0), NA_character_,
                               ifelse(LRae_0 >= lr, "Evidence for replication success",
                                      ifelse(LRae_0 <= 1/lr, "Evidence for replication failure",
                                             "Inconclusive evidence")))

    outcome <- ifelse(is.na(LRae_0), NA_character_,
                      ifelse(LRae_0 >= lr, "Replication success",
                             ifelse(LRae_0 <= 1/lr, "Replication failure",
                                    "Inconclusive")))

    data.frame(outcome = outcome, outcome_detailed = outcome_detailed,
               quantitative_value = LRae_0)
  }


  # Equivalence test ----
  assess_equivalence <- function(df, h = 0.2, alpha = 0.05){

    # h = 0.2 is defined in CORRELATION units (see assess_correspondence).
    # cor scale     -> es/se are already raw r, used as-is.
    # fisherz scale -> back-transform to r first so h keeps its meaning.
    needs_backtransform <- unique(df$scale) == "fisherz"

    compute_equiv <- function(es_o, se_o, es_r, se_r){
      if(any(is.na(c(es_o, se_o, es_r, se_r))) ||
         any(!is.finite(c(es_o, se_o, es_r, se_r))) ||
         se_o <= 0 || se_r <= 0){
        return(NA_real_)
      }
      if (needs_backtransform) {
        o <- fisherz_to_r(es_o, se_o)
        r <- fisherz_to_r(es_r, se_r)
        es_o <- o$r
        se_o <- o$se_r
        es_r <- r$r
        se_r <- r$se_r
      }
      tryCatch(equiv_test(es_o, es_r, se_o, se_r, h = h, alpha = alpha)$pval,
               error = function(e) NA_real_)
    }

    p_equiv <- mapply(compute_equiv, es_o = df$es_o, se_o = df$se_o,
                      es_r = df$es_r, se_r = df$se_r)

    # TOST: equivalent if both one-sided p <= alpha, i.e. max(p) <= alpha
    outcome_detailed <- ifelse(is.na(p_equiv), NA_character_,
                               ifelse(p_equiv <= alpha,
                                      "Original and replication effects are equivalent within threshold",
                                      "Original and replication effects are not demonstrably equivalent"))

    outcome <- ifelse(is.na(p_equiv), NA_character_,
                      ifelse(p_equiv <= alpha, "Replication success", "Replication failure"))

    data.frame(outcome = outcome, outcome_detailed = outcome_detailed,
               quantitative_value = p_equiv)
  }

  # Correspondence test ----
  assess_correspondence <- function(df, h = 0.2, alpha = 0.05){

    # The equivalence threshold h = 0.2 is defined in CORRELATION units
    # (a shift in corr. interpretation). The test must
    # therefore run on the correlation scale:
    #  * df$scale == "cor"     -> es/se are already raw r, use them as-is
    #  * df$scale == "fisherz" -> back-transform es/se to r (tanh + delta
    #    method) before testing, so h keeps its intended meaning.
    needs_backtransform <- unique(df$scale) == "fisherz"

    compute_corresp <- function(es_o, se_o, es_r, se_r){
      if(any(is.na(c(es_o, se_o, es_r, se_r))) ||
         any(!is.finite(c(es_o, se_o, es_r, se_r))) ||
         se_o <= 0 || se_r <= 0){
        return(list(ct = NA_character_, p_equiv = NA_real_))
      }
      if (needs_backtransform) {
        o <- fisherz_to_r(es_o, se_o)
        r <- fisherz_to_r(es_r, se_r)
        es_o <- o$r
        se_o <- o$se_r
        es_r <- r$r
        se_r <- r$se_r
      }
      tryCatch(corresp_test(es_o, es_r, se_o, se_r, h = h, alpha = alpha),
               error = function(e) list(ct = NA_character_, p_equiv = NA_real_))
    }

    res <- mapply(compute_corresp, es_o = df$es_o, se_o = df$se_o,
                  es_r = df$es_r, se_r = df$se_r, SIMPLIFY = FALSE)

    ct      <- vapply(res, function(x) x$ct,      character(1))
    p_equiv <- vapply(res, function(x) x$p_equiv, numeric(1))

    ### Success criterion
    ### Equivalence / Trivial Difference -> success; Difference -> failure;
    ### Indeterminacy -> NA (neither test conclusive)
    outcome <- ifelse(is.na(ct), NA_character_,
                      ifelse(ct %in% c("Equivalence", "Trivial Difference"), "Replication success",
                             ifelse(ct == "Difference", "Replication failure",
                                    ifelse(ct == "Indeterminacy", "Inconclusive", NA_character_))))

    data.frame(outcome = outcome, outcome_detailed = ct,
               quantitative_value = p_equiv)
  }

  # Minimum-effect test ----
  assess_minimum_effect <- function(df, h = NULL, shrinkage=0.25, alpha = 0.025){

    # Minimum-effect testing is the opposite (mirror) of equivalence testing (TOST):
    # H0 is rejected when the (1 - 2*alpha) CI for the difference
    # lies entirely outside [-h, h] --
    # The combined p-value is the minimum of the two one-sided p-values
    # while TOST takes the maximum.
    needs_backtransform <- unique(df$scale) == "fisherz"

    compute_min_effect <- function(es_o, se_o, es_r, se_r){
      if(any(is.na(c(es_o, se_o, es_r, se_r))) ||
         any(!is.finite(c(es_o, se_o, es_r, se_r))) ||
         se_o <= 0 || se_r <= 0){
        return(NA_real_)
      }
      if (needs_backtransform) {
        o <- fisherz_to_r(es_o, se_o); r <- fisherz_to_r(es_r, se_r)
        es_o <- o$r; se_o <- o$se_r; es_r <- r$r; se_r <- r$se_r
      }
      tryCatch({
        # orient both estimates to the original's direction
        flip  <- if (es_o < 0) -1 else 1
        es_o  <- es_o * flip
        es_r  <- es_r * flip

        h_use <- if (is.null(h)) (1 - shrinkage) * es_o else h
        if (is.na(h_use) || !is.finite(h_use)) return(NA_real_)

        # one-sided: evidence that theta_r lies BELOW h
        stats::pnorm((es_r - h_use) / se_r, lower.tail = TRUE)
      }, error = function(e) NA_real_)
    }

    p_min <- mapply(compute_min_effect, es_o = df$es_o, se_o = df$se_o,
                    es_r = df$es_r, se_r = df$se_r)

    rejected <- !is.na(p_min) & (p_min <= alpha)

    ### Success criterion: rejecting H0 = significant shrinkage below h = failure
    success <- !is.na(p_min) & !rejected

    outcome_detailed <- ifelse(is.na(p_min), NA_character_,
                               ifelse(rejected,
                                      "Replication effect is significantly smaller than the minimum effect of interest",
                                      "Replication effect is not significantly smaller than the minimum effect of interest"))

    outcome <- ifelse(is.na(p_min), NA_character_,
                      ifelse(success, "Replication success", "Replication failure"))

    data.frame(outcome = outcome, outcome_detailed = outcome_detailed,
               quantitative_value = p_min)
  }

  # Switching the names of the functions to metric argument
  metric_fns <- list(
    significance_r   = assess_rep_significance,
    significance_agg = assess_agg_significance,
    sum_p            = assess_sum_p,
    consistency_cir  = assess_or_consistency_with_rep,
    consistency_cio  = assess_rep_consistency_with_or,
    consistency_pi   = assess_rep_consistency_with_pi,
    small_telescopes = assess_small_telescopes,
    BF_JZS           = assess_BFJZS,
    BF_meta          = assess_BFmeta,
    BFe              = assess_BFe,
    BFr              = assess_BFr,
    BFs              = assess_BFs,
    pScept           = assess_pScept,
    snaphybrid       = assess_snaphybrid,
    Likratio         = assess_Likratio,
    equivalence      = assess_equivalence,
    correspondence   = assess_correspondence,
    minimum_effect   = assess_minimum_effect
  )

  unknown <- setdiff(metric, names(metric_fns))
  if (length(unknown) > 0) {
    stop("Metric(s) not recognised: ", paste(unknown, collapse=", "),
         "\nAvailable metrics: ", paste(names(metric_fns), collapse=", "))
  }

  # Metric-specific arguments (thresholds etc.) passed on to the metric function
  if (!is.null(metric_args)) {
    if (!is.list(metric_args) || is.null(names(metric_args)) ||
        any(!nzchar(names(metric_args))) || anyDuplicated(names(metric_args))) {
      stop("`metric_args` must be a named list with one element per metric, ",
           "e.g. list(BFe = list(gamma = 3)).")
    }
    not_requested <- setdiff(names(metric_args), metric)
    if (length(not_requested) > 0) {
      stop("`metric_args` given for metric(s) not in `metric`: ",
           paste(not_requested, collapse = ", "))
    }
    for (m in names(metric_args)) {
      a <- metric_args[[m]]
      if (!is.list(a) || (length(a) > 0 && (is.null(names(a)) || any(!nzchar(names(a)))))) {
        stop("`metric_args$", m, "` must be a named list of arguments.")
      }
      allowed <- setdiff(names(formals(metric_fns[[m]])), "df")
      unknown_args <- setdiff(names(a), allowed)
      if (length(unknown_args) > 0) {
        stop("Unknown argument(s) for metric '", m, "': ",
             paste(unknown_args, collapse = ", "),
             ". Available: ", paste(allowed, collapse = ", "))
      }
    }
  }

  results <- lapply(metric, function(m) {
    res <- do.call(metric_fns[[m]], c(list(df), metric_args[[m]]))
    res$study_id <- df$study_id
    res$metric   <- m
    res
  })

  long <- do.call(rbind, results)
  long <- long[order(long$study_id, long$metric), ]
  rownames(long) <- NULL

  wide <- reshape(
    long[, c("study_id", "metric", "quantitative_value", "outcome", "outcome_detailed")],
    idvar     = "study_id",
    timevar   = "metric",
    direction = "wide"
  )
  rownames(wide) <- NULL

  names(wide) <- gsub("^quantitative_value\\.(.+)$","\\1_value",   names(wide))
  names(wide) <- gsub("^outcome\\.(.+)$",           "\\1_outcome", names(wide))
  names(wide) <- gsub("^outcome_detailed\\.(.+)$",  "\\1_outcome_detailed", names(wide))

  bf_metrics <- c("BF_JZS", "BF_meta", "BFe", "BFr", "BFs", "snaphybrid")
  for (m in metric) {
    wide[[paste0(m, "_type")]] <- ifelse(m %in% bf_metrics, "Bayes factor", "p-value")
  }

  metric_cols <- unlist(lapply(metric, function(m) {
    c(paste0(m, "_type"), paste0(m, "_value"), paste0(m, "_outcome"), paste0(m, "_outcome_detailed"))
  }))
  wide <- wide[, c("study_id", metric_cols)]

  data_df <- df[, c("study_id", "es_o", "se_o", "es_r", "se_r")]
  rownames(data_df) <- NULL

  return(list(data=data_df, results=wide))
}

# Function 3: Get replication results -----------------------------------------

#' get_replication_results: Get replication results
#'
#' Extracts the qualitative outcomes, quantitative values, or prepared
#' effect-size data from an object returned by
#' [assess_replication_outcome()].
#'
#' @param x A replication-assessment object returned by
#'   [assess_replication_outcome()]. It must be a list containing components
#'   named `data` and `results`.
#' @param type Character string containing which information to extract from the results.
#'   Options are:
#'   \describe{
#'     \item{`"qualitative"`}{Returns the overall and detailed replication
#'       outcomes for each selected metric. This is the default.}
#'     \item{`"quantitative"`}{Returns the quantitative value and value-type
#'       label for each selected metric.}
#'     \item{`"data"`}{Returns the prepared original and replication effect
#'       sizes used in the assessment.}
#'   }
#'
#' @details
#' This function does not recalculate any replication metric. It provides a
#' convenient way to extract different parts of the output produced by
#' [assess_replication_outcome()].
#'
#' For qualitative results, the columns are arranged by metric as the overall
#' outcome followed by its detailed interpretation. For quantitative results,
#' the columns are arranged as the quantitative value followed by its
#' corresponding type.
#'
#' @return A data frame with one row per original-replication pair:
#'   \describe{
#'     \item{If `type = "qualitative"`}{Returns `study_id` followed by
#'       `<metric>_outcome` and `<metric>_outcome_detailed` for every assessed
#'       metric.}
#'     \item{If `type = "quantitative"`}{Returns `study_id` followed by
#'       `<metric>_value` and `<metric>_type` for every assessed metric.}
#'     \item{If `type = "data"`}{Returns `study_id`, `es_o`, `se_o`, `es_r`,
#'       and `se_r`, representing the effect estimates and standard errors
#'       used in the assessment.}
#'   }
#'
#' @seealso [prepare_effect_sizes()], [assess_replication_outcome()]
#'
#' @export
#'
#' @examples
#' df <- prepare_effect_sizes(
#'   r_o = c(0.30, 0.20),
#'   n_o = c(100, 150),
#'   r_r = c(0.25, 0.05),
#'   n_r = c(120, 180),
#'   effectsize = "COR",
#'   scale = "fisherz"
#' )
#'
#' replication_results <- assess_replication_outcome(df,
#'        metric = c("significance_r",
#'                   "consistency_pi",
#'                   "BF_JZS"))
#'
#' # Extract the replication classifications
#' get_replication_results(
#'   replication_results,
#'   type = "qualitative"
#' )
#'
#' # Extract the metric values
#' get_replication_results(
#'   replication_results,
#'   type = "quantitative"
#' )
#'
#' # Extract the effect sizes used in the analysis
#' get_replication_results(
#'   replication_results,
#'   type = "data"
#' )
get_replication_results <- function(x, type = c("qualitative", "quantitative", "data")) {

  type <- match.arg(type)

  stopifnot(
    is.list(x),
    all(c("data", "results") %in% names(x))
  )

  if (type == "data") {
    return(x$data)
  }

  if (type == "qualitative") {
    outcome_cols          <- grep("_outcome$",          names(x$results), value = TRUE)
    outcome_detailed_cols <- grep("_outcome_detailed$", names(x$results), value = TRUE)
    # Interleave: metric1_outcome, metric1_outcome_detailed, metric2_outcome, ...
    interleaved <- c(rbind(outcome_cols, outcome_detailed_cols))
    return(x$results[, c("study_id", interleaved)])
  }

  if (type == "quantitative") {
    value_cols <- grep("_value$", names(x$results), value = TRUE)
    type_cols  <- grep("_type$",  names(x$results), value = TRUE)
    interleaved <- c(rbind(value_cols, type_cols))
    return(x$results[, c("study_id", interleaved)])
  }
}
