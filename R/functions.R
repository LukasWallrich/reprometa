
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
#' Standardized mean differences are first converted to correlations and then
#' to Fisher-\eqn{z} correlations. Log odds ratios are converted successively
#' to standardized mean differences, correlations, and Fisher-\eqn{z}
#' correlations. Inputs supplied as `"ZCOR"` are retained on their supplied
#' scale.
#'
#' Optional reported confidence intervals, p-values, and pattern indicators
#' are attached to the resulting data frame.
#'
#' @return A data frame with one row per original–replication pair. It contains
#'   at least the following columns:
#'   \describe{
#'     \item{`study_id`}{Sequential identifier for the study pair.}
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
#' # Retain correlations on their original scale
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
                                 # --- Optional arguments: reported p-values + coding metadata.
                                 #     orig_p_value is kept only when orig_p_value_type =="exact" (if type is supplied),
                                 #     doubled when reported one-tailed (if tails is supplied)
                                 #     with limits to (0,1] if correlation.
                                 #     These provided p-values are used by sum_p and by pScept via the p2z() route.
                                 orig_p_value=NULL, orig_p_value_type=NULL,
                                 orig_p_value_tails=NULL, repli_p_value=NULL,
                                 # --- Optional: test-statistic metadata for exact, test-type-dependent SE.
                                 #     When supplied, se_o/se_r are overwritten: on the cor scale with
                                 #     sqrt((1-r^2)/df1) for t-type (else /size); on the
                                 #     fisherz scale with 1/sqrt(n_eff - 3), where n_eff is
                                 #     the same effective df/size.
                                 stat_type_o=NULL, dof_o=NULL,
                                 effective_df_o=NULL, effective_sample_size_o=NULL,
                                 stat_type_r=NULL, dof_r=NULL,
                                 effective_df_r=NULL, effective_sample_size_r=NULL,
                                 # --- Optional: qualitative direction/pattern flag.
                                 #     Used to gate pScept / BFr on the cor scale
                                 #     Without this argument these metrics fall back to
                                 #     a sign-match approximation.
                                 pattern_met=NULL) {

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

  # Same effective df/size concept translated to the Fisher-z SE, 1/sqrt(n_eff-3)
  compute_se_score_z <- function(r, stat_type, dof, effective_df,
                                 sample_size, effective_sample_size) {
    if (is.na(r) || is.na(stat_type)) return(NA_real_)
    df1  <- ifelse(is.na(effective_df), dof, effective_df)
    size <- ifelse(is.na(effective_sample_size), sample_size, effective_sample_size)
    n_eff <- if (stat_type == "t" && !is.na(df1)) df1 else size
    if (is.na(n_eff) || n_eff <= 3) return(NA_real_)
    1 / sqrt(n_eff - 3)
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

      if (scale == "cor") {
        es_o_i <- r_o[i]
        es_r_i <- r_r[i]
        if (have_se_meta) {
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
        # Exact Fisher-z transform (no Hedges-Olkin bias correction on the standard
        # error): containment checks against atanh()-transformed CI bounds only if es is atanh(r) exactly.
        es_o_i <- atanh(r_o[i])
        es_r_i <- atanh(r_r[i])
        if (have_se_meta) {
          # effective df/size concept translated to the Fisher-z scale
          se_o_i <- compute_se_score_z(r_o[i], pad(stat_type_o)[i], pad(dof_o)[i],
                                       pad(effective_df_o)[i], n_o[i],
                                       pad(effective_sample_size_o)[i])
          se_r_i <- compute_se_score_z(r_r[i], pad(stat_type_r)[i], pad(dof_r)[i],
                                       pad(effective_df_r)[i], n_r[i],
                                       pad(effective_sample_size_r)[i])
        } else {
          se_o_i <- sqrt(1 / n_o[i] + (6 - r_o[i]^2) / (2 * n_o[i]^2))
          se_r_i <- sqrt(1 / n_r[i] + (6 - r_r[i]^2) / (2 * n_r[i]^2))
        }
      }

      data.frame(es_o=es_o_i, se_o=se_o_i, es_r=es_r_i, se_r=se_r_i,
                 n_o=n_o[i], n_r=n_r[i],
                 n1_o=NA_real_, n2_o=NA_real_, n1_r=NA_real_, n2_r=NA_real_)

    } else if (es == "SMD") {
      if (is.null(d_o) || is.null(d_r) || is.null(se_o) || is.null(se_r) ||
          is.null(n1_o) || is.null(n2_o) || is.null(n1_r) || is.null(n2_r))
        stop("effectsize 'SMD' requires d_o, d_r, se_o, se_r, n1_o, n2_o, n1_r, n2_r - missing for study ", i)

      es_o_i <- r_to_zr(d_to_r(d_o[i], se_o[i], n1=n1_o[i], n2=n2_o[i])[1],
                        d_to_r(d_o[i], se_o[i], n1=n1_o[i], n2=n2_o[i])[2])[["Fisher Z correlation"]]
      es_r_i <- r_to_zr(d_to_r(d_r[i], se_r[i], n1=n1_r[i], n2=n2_r[i])[1],
                        d_to_r(d_r[i], se_r[i], n1=n1_r[i], n2=n2_r[i])[2])[["Fisher Z correlation"]]
      se_o_i <- r_to_zr(d_to_r(d_o[i], se_o[i], n1=n1_o[i], n2=n2_o[i])[1],
                        d_to_r(d_o[i], se_o[i], n1=n1_o[i], n2=n2_o[i])[2])[["Std. Err."]]
      se_r_i <- r_to_zr(d_to_r(d_r[i], se_r[i], n1=n1_r[i], n2=n2_r[i])[1],
                        d_to_r(d_r[i], se_r[i], n1=n1_r[i], n2=n2_r[i])[2])[["Std. Err."]]

      data.frame(es_o=es_o_i, se_o=se_o_i, es_r=es_r_i, se_r=se_r_i,
                 n_o=n1_o[i]+n2_o[i], n_r=n1_r[i]+n2_r[i],
                 n1_o=n1_o[i], n2_o=n2_o[i], n1_r=n1_r[i], n2_r=n2_r[i])

    } else if (es == "LOGOR") {
      if (is.null(logor_o) || is.null(logor_r) || is.null(se_o) || is.null(se_r) ||
          is.null(n1_o) || is.null(n2_o) || is.null(n1_r) || is.null(n2_r))
        stop("effectsize 'LOGOR' requires logor_o, logor_r, se_o, se_r, n1_o, n2_o, n1_r, n2_r - missing for study ", i)

      d_o_i  <- logor_to_d(logor=logor_o[i], se_logor=se_o[i])[["Cohen's d"]]
      d_r_i  <- logor_to_d(logor=logor_r[i], se_logor=se_r[i])[["Cohen's d"]]
      se_o_i <- logor_to_d(logor=logor_o[i], se_logor=se_o[i])[["Std. Err."]]
      se_r_i <- logor_to_d(logor=logor_r[i], se_logor=se_r[i])[["Std. Err."]]

      es_o_i <- r_to_zr(d_to_r(d_o_i, se_o_i, n1=n1_o[i], n2=n2_o[i])[1],
                        d_to_r(d_o_i, se_o_i, n1=n1_o[i], n2=n2_o[i])[2])[["Fisher Z correlation"]]
      es_r_i <- r_to_zr(d_to_r(d_r_i, se_r_i, n1=n1_r[i], n2=n2_r[i])[1],
                        d_to_r(d_r_i, se_r_i, n1=n1_r[i], n2=n2_r[i])[2])[["Fisher Z correlation"]]
      se_o_i <- r_to_zr(d_to_r(d_o_i, se_o_i, n1=n1_o[i], n2=n2_o[i])[1],
                        d_to_r(d_o_i, se_o_i, n1=n1_o[i], n2=n2_o[i])[2])[["Std. Err."]]
      se_r_i <- r_to_zr(d_to_r(d_r_i, se_r_i, n1=n1_r[i], n2=n2_r[i])[1],
                        d_to_r(d_r_i, se_r_i, n1=n1_r[i], n2=n2_r[i])[2])[["Std. Err."]]

      data.frame(es_o=es_o_i, se_o=se_o_i, es_r=es_r_i, se_r=se_r_i,
                 n_o=n1_o[i]+n2_o[i], n_r=n1_r[i]+n2_r[i],
                 n1_o=n1_o[i], n2_o=n2_o[i], n1_r=n1_r[i], n2_r=n2_r[i])

    } else if (es == "ZCOR") {
      if (is.null(z_o) || is.null(z_r) || is.null(se_o) || is.null(se_r) || is.null(n_o) || is.null(n_r))
        stop("effectsize 'ZCOR' requires z_o, z_r, se_o, se_r, n_o, n_r - missing for study ", i)

      data.frame(es_o = z_o[i], se_o = se_o[i], es_r = z_r[i], se_r = se_r[i],
                 n_o = n_o[i], n_r = n_r[i],
                 n1_o = NA_real_, n2_o = NA_real_, n1_r = NA_real_, n2_r = NA_real_)
    }
  }

  df <- do.call(rbind, lapply(seq_len(n_studies), convert_one))
  df <- cbind(study_id=seq_len(n_studies), effectsize=effectsize, scale=scale, df)
  rownames(df) <- NULL
  df <- df[, colSums(is.na(df)) < nrow(df)]

  # --- Attach precomputed CI bounds on the analysis scale ---------------------
  # Input bounds are on the correlation scale
  # scale="cor"    -> attached raw
  # scale="fisherz"-> atanh()-transformed (bounds clamped just inside (-1,1) so
  #                   atanh stays finite), under es_lb_*/es_ub_* names.
  # atanh is strictly monotonic, so CI-based assessments against the
  # transformed bounds are identical to raw-r decisions.
  if (!is.null(orig_conv_r_lb) || !is.null(repli_conv_r_lb)) {
    z_bound <- function(b) atanh(pmax(pmin(b, 1 - 1e-12), -1 + 1e-12))
    if (scale == "cor") {
      if (!is.null(orig_conv_r_lb))  df$orig_conv_r_lb  <- orig_conv_r_lb
      if (!is.null(orig_conv_r_ub))  df$orig_conv_r_ub  <- orig_conv_r_ub
      if (!is.null(repli_conv_r_lb)) df$repli_conv_r_lb <- repli_conv_r_lb
      if (!is.null(repli_conv_r_ub)) df$repli_conv_r_ub <- repli_conv_r_ub
    } else {
      if (!is.null(orig_conv_r_lb))  df$es_lb_o <- z_bound(orig_conv_r_lb)
      if (!is.null(orig_conv_r_ub))  df$es_ub_o <- z_bound(orig_conv_r_ub)
      if (!is.null(repli_conv_r_lb)) df$es_lb_r <- z_bound(repli_conv_r_lb)
      if (!is.null(repli_conv_r_ub)) df$es_ub_r <- z_bound(repli_conv_r_ub)
    }
  }

  # --- Attach standardized p-values (scale-free) ------------------------------
  #  only "exact"-type original p-values are usable (if type metadata supplied)
  #  one-tailed original p-values are doubled to two-tailed (if tails supplied)
  #  clamped to (0,1] so ReplicationSuccess::p2z() accepts them
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

  # --- Attach qualitative direction/pattern flag ------------------------------
  if (!is.null(pattern_met)) {
    df$repli_pattern_criteria_met <- pattern_met
  }

  return(df)
}

# Function 2: Assess replication success criteria ----------------------------

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
#'       criterion.}
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
#'   overall type-I error rate of
#'   \eqn{\alpha^2 = 0.025^2 = 0.000625} for metrics for which a calibrated
#'   version is implemented. Specifically:
#'   \itemize{
#'     \item `"significance_agg"` uses an unweighted Stouffer combination with
#'       a threshold calibrated to \eqn{\alpha^2};
#'     \item `"BF_JZS"` uses a Bayes-factor threshold of 2.35;
#'     \item `"BFs"` uses a Bayes-factor threshold of 4.08; and
#'     \item `"pScept"` uses the controlled rather than the golden
#'       recalibration.
#'   }
#'   Metrics already calibrated by construction, including the two-trials rule
#'   and sum-of-p-values criterion, are unaffected. Other metrics retain their
#'   conventional thresholds.
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
assess_replication_outcome <- function(df, metric, calibrated=TRUE) {

  # `calibrated = TRUE` switches the metrics whose conventional thresholds do
  # NOT control the overall type-I error rate at alpha^2 = 0.025^2 = 0.000625
  # over to calibrated versions:
  #   BF_JZS  : gamma 1 -> 2.35
  #   BFs     : gamma 1 -> 4.08
  #   pScept  : golden -> controlled recalibration
  #   meta    : fixed-effect pooling -> unweighted Stouffer at sqrt(2)*z_{1-alpha^2}
  # Metrics already calibrated at alpha^2 by construction (two-trials rule,
  # sum of p-values) are unaffected.

  stopifnot(is.logical(calibrated), length(calibrated) == 1)

  stopifnot(is.data.frame(df), is.character(metric))
  stopifnot(all(c("es_o", "se_o", "es_r", "se_r", "study_id") %in% names(df)))
  stopifnot(metric %in% c("significance_r","significance_agg","sum_p","consistency_cir",
                          "consistency_cio","consistency_pi","small_telescopes",
                          "BF_JZS","BFe","BFr","BFs","pScept","snaphybrid",
                          "equivalence","correspondence","Likratio","minimum_effect"))


  # Significance criterion (vote counting) ----
  assess_rep_significance <- function(df, alpha = 0.025) {

    # Sidedness follows the analysis scale:
    #  * fisherz -> ONE-SIDED p in the direction of the estimate, compared to
    #    alpha (original implementation).
    #  * cor     -> TWO-SIDED p compared to 2*alpha (numerically the same
    #    decision when combined with the same-direction requirement below)
    if (unique(df$scale) == "cor") {
      p_o <- 2 * pnorm(abs(df$es_o / df$se_o), lower.tail = FALSE)
      p_r <- 2 * pnorm(abs(df$es_r / df$se_r), lower.tail = FALSE)
      alpha <- 2 * alpha
    } else {
      p_o <- zr_to_p(zr = df$es_o, se = df$se_o)
      p_r <- zr_to_p(zr = df$es_r, se = df$se_r)
    }

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
  assess_agg_significance <- function(df, alpha = 0.025) {

    if (calibrated) {
      # --- CALIBRATED: unweighted Stouffer, exact overall T1E at alpha^2 -----

      z_o <- df$es_o / df$se_o
      z_r <- df$es_r / df$se_r
      flip <- ifelse(!is.na(z_o) & z_o < 0, -1, 1)   # orient to the original's direction
      z_o <- z_o * flip
      z_r <- z_r * flip

      z_sum  <- z_o + z_r
      thresh <- sqrt(2) * stats::qnorm(1 - alpha^2)
      p_meta <- stats::pnorm(z_sum / sqrt(2), lower.tail = FALSE)   # combined one-sided p

      success <- !is.na(z_sum) & is.finite(z_sum) & (z_sum >= thresh)

      outcome_detailed <- ifelse(is.na(z_sum) | !is.finite(z_sum), NA_character_,
                                 ifelse(success,
                                        "Aggregated effect exceeds the success threshold",
                                        "Aggregated effect does not exceed the success threshold"))
      outcome <- ifelse(is.na(z_sum) | !is.finite(z_sum), NA_character_,
                        ifelse(success, "Replication success", "Replication failure"))

      return(data.frame(outcome = outcome, outcome_detailed = outcome_detailed,
                        quantitative_value = p_meta))
    }


    # Calculates the pooled effect estimate and CI lower/upper bound
    effect_comparison <- compare_effectsizes(df$es_o, df$se_o, df$es_r, df$se_r)

    ### Calculate p-value for this metric
    p_meta <- effect_comparison$pval

    # Success criterion, matching repli_binary.R's meta_success exactly:
    #   sig = sign(lb_b) == sign(ub_b)      (pooled CI excludes zero)
    #   dir = sign(orig_conv_r) == sign(lb_b)  (direction matches the original)
    # NOTE: this needs signed es_o/ci_ub -- if you're feeding in abs()'d
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

  assess_sum_p <- function(df, alpha = 0.025, threshold = 0.025*sqrt(2)){

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

    both_p <- sum_p < threshold

    # both_p <- ifelse(sum_p <= 1,
    #                  sum_p^2 / 2 < alpha^2,
    #                  1 - (2 - sum_p)^2 / 2 < alpha^2)

    ### Success criterion: the sum-of-p test is direction-agnostic, so AND with
    ### the qualitative pattern flag -- exactly like repli_binary.R's
    ### `sum_p <- both_p & repli_pattern_criteria_met`
    pattern_met <- if ("repli_pattern_criteria_met" %in% names(df)) {
      df$repli_pattern_criteria_met
    } else {
      sign(df$es_o) == sign(df$es_r)
    }

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
  assess_or_consistency_with_rep <- function(df) {

    if (unique(df$scale) == "cor" && all(c("repli_conv_r_lb","repli_conv_r_ub") %in% names(df))) {

      outcome <- ifelse(df$es_o > df$repli_conv_r_lb & df$es_o < df$repli_conv_r_ub,
                        "Replication success", "Replication failure")

      data.frame(outcome=outcome, outcome_detailed=outcome, quantitative_value=NA_real_)

    } else if (unique(df$scale) == "fisherz" && all(c("es_lb_r","es_ub_r") %in% names(df))) {
      # Containment against supplied bounds on the z scale (atanh of the study's
      # reported CI). atanh is monotonic, so decision-identical to raw-r.
      outcome <- ifelse(is.na(df$es_o) | is.na(df$es_lb_r) | is.na(df$es_ub_r), NA_character_,
                  ifelse(df$es_o > df$es_lb_r & df$es_o < df$es_ub_r,
                         "Replication success", "Replication failure"))
      data.frame(outcome=outcome, outcome_detailed=outcome, quantitative_value=NA_real_)

    } else {

      p_cir <- pnorm((df$es_r - df$es_o) / df$se_r, mean = 0, sd = 1, lower.tail = FALSE)
      cis <- compute_ci(df$es_r, df$se_r)

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
  assess_rep_consistency_with_or <- function(df) {

    if (unique(df$scale) == "cor" && all(c("orig_conv_r_lb","orig_conv_r_ub") %in% names(df))) {
      # --- precomputed bounds directly
      outcome <- ifelse(df$es_r > df$orig_conv_r_lb & df$es_r < df$orig_conv_r_ub,
                        "Replication success", "Replication failure")
      data.frame(outcome=outcome, outcome_detailed=outcome, quantitative_value=NA_real_)

    } else if (unique(df$scale) == "fisherz" && all(c("es_lb_o","es_ub_o") %in% names(df))) {
      # Same monotonicity argument: decision-identical to raw-r.
      outcome <- ifelse(is.na(df$es_r) | is.na(df$es_lb_o) | is.na(df$es_ub_o), NA_character_,
                  ifelse(df$es_r > df$es_lb_o & df$es_r < df$es_ub_o,
                         "Replication success", "Replication failure"))
      data.frame(outcome=outcome, outcome_detailed=outcome, quantitative_value=NA_real_)

    } else {
      # --- existing Fisher-z / self-computed-CI path, unchanged ---
      p_cio <- pnorm((df$es_o - df$es_r) / df$se_o, mean = 0, sd = 1, lower.tail = FALSE)
      cis <- compute_ci(df$es_o, df$se_o)
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
  assess_rep_consistency_with_pi <- function(df) {

    ### Calculate p-value for this metric
    p_PI <- pnorm((df$es_o - df$es_r) / sqrt(df$se_o^2 + df$se_r^2), mean = 0, sd = 1, lower.tail = FALSE)

    ### Calculate prediction interval
    pred_intervals <- compute_prediction_interval(df$es_o, df$se_o, df$se_r)

    # On the correlation scale, cap the bounds to [-1, 1] -- matches
    # repli_binary.R exactly ("these can go above 1 and below -1, so just
    # manually capping them for now"). Not applicable on the Fisher-z scale,
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

  assess_small_telescopes <- function(df, power=1/3, alpha=0.025) {

    n1_o <- if ("n1_o" %in% names(df)) df$n1_o else rep(NA_real_, nrow(df))
    n2_o <- if ("n2_o" %in% names(df)) df$n2_o else rep(NA_real_, nrow(df))

    if (unique(df$scale) == "cor" && all(df$effectsize == "COR")) {
      small_effect <- mapply(
        function(r, n) {

          if (is.na(r) || !is.finite(r) || is.na(n) || !is.finite(n) || n < 4) {
            return(NA_real_)
          }

          thr <- tryCatch(
            pwr::pwr.r.test(n = n, power = power)$r,
            error = function(e) NA_real_
          )

          if (is.na(thr) || !is.finite(thr)) {
            return(NA_real_)
          }

          if (sign(r) == 1) {
            thr
          } else {
            -thr
          }
        }, r = df$es_o, n = df$n_o)

      # Compare against the precomputed replication CI bounds directly --
      # never recompute this CI, matching repli_binary.R exactly.
      outcome_detailed <- ifelse(sign(small_effect) == 1,
                                 ifelse(df$repli_conv_r_ub > small_effect,
                                        "original study had >= 33% power to detect replication effect",
                                        "original study had < 33% power to detect replication effect"),
                                 ifelse(small_effect > df$repli_conv_r_lb,
                                        "original study had >= 33% power to detect replication effect",
                                        "original study had < 33% power to detect replication effect"))
      quant_val <- small_effect

    } else if (unique(df$scale) == "fisherz" &&
               all(c("es_lb_r","es_ub_r") %in% names(df)) &&
               all(df$effectsize %in% c("COR","ZCOR"))) {
      # Path on the Fisher-z scale ---
      # threshold r33 from pwr.r.test with ITS OWN defaults (sig.level=0.05,
      # two.sided. atanh() transformed onto the z scale, compared against the supplied z-scale CI bounds.
      small_effect <- mapply(function(z, n) {
        if (is.na(z) || is.na(n) || !is.finite(n) || n < 4) return(NA_real_)
        thr <- tryCatch(atanh(pwr::pwr.r.test(n = n, power = power)$r),
                        error = function(e) NA_real_)
        if (is.na(thr) || !is.finite(thr)) return(NA_real_)
        if (sign(z) == 1) thr else -thr
      }, z = df$es_o, n = df$n_o)

      hit <- ifelse(sign(small_effect) == 1,
                    df$es_ub_r > small_effect,
                    small_effect > df$es_lb_r)
      outcome_detailed <- ifelse(is.na(hit), NA_character_,
                          ifelse(hit,
                                 "original study had >= 33% power to detect replication effect",
                                 "original study had < 33% power to detect replication effect"))
      quant_val <- small_effect

    } else {
      # --- existing Fisher-z / probabilistic path, unchanged ---
      es_with_33power <- mapply(function(effectsize, n1, n2, n) {
        if (effectsize %in% c("SMD", "LOGOR")) {
          if (is.na(n1) || is.na(n2)) return(NA_real_)
          d33 <- find_z_for_power(n1 = n1, n2 = n2, power = power, alpha = alpha)
          r_to_zr(d_to_r(d33, n1 = n1, n2 = n2)[["r"]])[[1]]
        } else {
          if (is.na(n) || n < 4) return(NA_real_)
          power.z.test(power = power, n = n-3, type = "one.sample",
                       sig.level = alpha, alternative = "one.sided")$delta
        }
      }, effectsize=df$effectsize, n1=n1_o, n2=n2_o, n=df$n_o)
      es_with_33power <- as.numeric(es_with_33power)

      p_telescopes <- pnorm((df$es_r - es_with_33power) / df$se_r, mean=0, sd=1, lower.tail=TRUE)
      outcome_detailed <- ifelse(p_telescopes >= alpha,
                                 "original study had >= 33% power to detect replication effect",
                                 "original study had < 33% power to detect replication effect")
      quant_val <- p_telescopes
    }

    outcome <- ifelse(outcome_detailed == "original study had >= 33% power to detect replication effect",
                      "Replication success", "Replication failure")

    data.frame(outcome=outcome, outcome_detailed=outcome_detailed, quantitative_value = quant_val)
  }

  # JZS Bayes Factor ----
  assess_BFJZS <- function(df, gamma = 1){

    # Calibrated threshold: gamma = 2.35 gives overall T1E = alpha^2
    if (calibrated) gamma <- 2.35

    # Computes integrand of the likelihood and the Cauchy prior (JZS prior)
    integrate_jzs <- function(theta, es, se){
      stats::dnorm(x = es, mean = theta, sd = se) *
        stats::dcauchy(theta, location = 0, scale = 1)
    }

    # Function for computing vectorized integral. Divide by marginal likelihood under the null to obtain BF
    compute_bfjzs <- function(es, se) {

      # Handle missing or invalid inputs
      if (is.na(es) || is.na(se) ||
          !is.finite(es) || !is.finite(se) ||
          se <= 0) {
        return(NA_real_)
      }

      # Numerical integration with error handling
      bf <- tryCatch({

        log(stats::integrate(integrate_jzs, lower = -Inf, upper = Inf,
                             es = es, se = se)$value) -
          log(stats::dnorm(x = es, mean = 0, sd = se))
      }, error = function(e) {
        NA_real_
      })

      bf
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

    # --- gate with the real direction/pattern flag, matching repli_binary.R's
    #     `bayes_rep <- repli_pattern_criteria_met & r_bf < 1` ---
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

    # Calibrated threshold: gamma = 4.08 gives overall T1E = alpha^2
    if (calibrated) gamma <- 4.08

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

  assess_pScept <- function(df, alpha = 0.025){

    # Sidedness follows the analysis scale:
    #  * scale == "cor"     -> TWO-SIDED: z's from the reported
    #    p-values via p2z(), threshold levelSceptical(0.05, "two.sided",
    #    "golden") = 0.1234 (hardcoded 0.05 to match repli_binary.R, which
    #    calls it that way independent of any general alpha), success gated
    #    by the qualitative pattern flag (two.sided pSceptical doesn't check
    #    direction itself).
    #  * scale == "fisherz" -> ONE-SIDED, original implementation: z's from
    #    the effect size / SE ratio, threshold levelSceptical(alpha,
    #    "one.sided", "golden"); no pattern gate (one.sided pSceptical
    #    handles direction internally).
    score_mode <- unique(df$scale) == "cor"

    # Calibrated: swap the "golden" recalibration for the "controlled" one,
    # which achieves EXACT overall T1E control at alpha^2 by making the success
    # threshold a function of the variance ratio c = n_r/n_o (so unlike golden,
    # the threshold varies per study pair and must be computed row-wise).
    ps_type <- if (calibrated) "controlled" else "golden"

    # controlled levelSceptical() requires a finite c for every element, so
    # compute it only where c is usable and leave NA elsewhere (those rows are
    # already NA in pScept and drop out of the success test regardless).
    controlled_level <- function(cc, alt) {
      out <- rep(NA_real_, length(cc))
      ok  <- !is.na(cc) & is.finite(cc) & cc > 0
      if (any(ok)) {
        out[ok] <- ReplicationSuccess::levelSceptical(level = alpha, c = cc[ok],
                                                      alternative = alt, type = "controlled")
      }
      out
    }


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

    # --- guard against dropped columns (COR-only data has no n1_o/n2_o/n1_r/n2_r) ---
    n1_o <- if ("n1_o" %in% names(df)) df$n1_o else rep(NA_real_, nrow(df))
    n2_o <- if ("n2_o" %in% names(df)) df$n2_o else rep(NA_real_, nrow(df))
    n1_r <- if ("n1_r" %in% names(df)) df$n1_r else rep(NA_real_, nrow(df))
    n2_r <- if ("n2_r" %in% names(df)) df$n2_r else rep(NA_real_, nrow(df))

    # --- resolve n_o/n_r ONCE, vectorized ---
    n_o_resolved <- ifelse(!is.na(n1_o) & !is.na(n2_o), n1_o + n2_o, df$n_o)
    n_r_resolved <- ifelse(!is.na(n1_r) & !is.na(n2_r), n1_r + n2_r, df$n_r)

    if (score_mode) {
      p_o <- if ("orig_p_value"  %in% names(df)) df$orig_p_value  else rep(NA_real_, nrow(df))
      p_r <- if ("repli_p_value" %in% names(df)) df$repli_p_value else rep(NA_real_, nrow(df))

      pScept <- mapply(compute_pScept_cor,
                       es_o = df$es_o, se_o = df$se_o, es_r = df$es_r, se_r = df$se_r,
                       n_o = n_o_resolved, n_r = n_r_resolved,
                       p_o = p_o, p_r = p_r)

      # golden: one fixed threshold. controlled: threshold depends on c = n_r/n_o,
      # so it is computed per study pair (levelSceptical is vectorised over c).
      skep_level <- if (calibrated) {
        controlled_level(n_r_resolved / n_o_resolved, alt = "one.sided")
      } else {
        ReplicationSuccess::levelSceptical(level = 0.05,
                                           alternative = "two.sided",
                                           type = "golden")
      }

      pattern_met <- if ("repli_pattern_criteria_met" %in% names(df)) {
        df$repli_pattern_criteria_met
      } else {
        sign(df$es_o) == sign(df$es_r)
      }

      success <- !is.na(pScept) & !is.na(pattern_met) & !is.na(skep_level) & (pScept < skep_level) & pattern_met
      na_mask <- is.na(pScept) | is.na(pattern_met)

    } else {
      pScept <- mapply(compute_pScept_z,
                       es_o = df$es_o, se_o = df$se_o, es_r = df$es_r, se_r = df$se_r,
                       n_o = n_o_resolved, n_r = n_r_resolved)

      skep_level <- if (calibrated) {
        controlled_level(n_r_resolved / n_o_resolved, "one.sided")
      } else {
        ReplicationSuccess::levelSceptical(level = alpha,
                                           alternative = "one.sided",
                                           type = "golden")
      }

      success <- !is.na(pScept) & !is.na(skep_level) & (pScept <= skep_level)
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

    # compute_pScept_cor <- function(es_r, se_r, es_o, se_o, n_o, n_r,
    #                                p_o = NA_real_, p_r = NA_real_){
    #   if (any(is.na(c(n_o, n_r))) || n_o <= 0 || n_r <= 0) return(NA_real_)
    #   p_ok <- function(p) !is.na(p) && is.finite(p) && p > 0 && p <= 1
    #   if (p_ok(p_o) && p_ok(p_r)) {
    #     zo <- ReplicationSuccess::p2z(p_o)
    #     zr <- ReplicationSuccess::p2z(p_r)
    #   } else {
    #     return(NA_real_)   # no valid exact p -> unusable -> NA
    #   }
    #   tryCatch(
    #     ReplicationSuccess::pSceptical(zo = zo, zr = zr, c = n_r/n_o,
    #                                    alternative = "two.sided", type = "golden"),
    #     error = function(e) NA_real_
    #   )
    # }
    #
    # compute_pScept_z <- function(es_r, se_r, es_o, se_o, n_o, n_r){
    #   if (any(is.na(c(es_r, se_r, es_o, se_o, n_o, n_r))) ||
    #       any(!is.finite(c(es_r, se_r, es_o, se_o, n_o, n_r))) ||
    #       se_o <= 0 || se_r <= 0 || n_o <= 0 || n_r <= 0){
    #     return(NA_real_)
    #   }
    #   tryCatch(
    #     ReplicationSuccess::pSceptical(zo = es_o/se_o, zr = es_r/se_r, c = n_r/n_o,
    #                                    alternative = "one.sided", type = "golden"),
    #     error = function(e) NA_real_
    #   )
    # }
    #
    # # --- guard against dropped columns (COR-only data has no n1_o/n2_o/n1_r/n2_r) ---
    # n1_o <- if ("n1_o" %in% names(df)) df$n1_o else rep(NA_real_, nrow(df))
    # n2_o <- if ("n2_o" %in% names(df)) df$n2_o else rep(NA_real_, nrow(df))
    # n1_r <- if ("n1_r" %in% names(df)) df$n1_r else rep(NA_real_, nrow(df))
    # n2_r <- if ("n2_r" %in% names(df)) df$n2_r else rep(NA_real_, nrow(df))
    #
    # # --- resolve n_o/n_r ONCE, vectorized ---
    # n_o_resolved <- ifelse(!is.na(n1_o) & !is.na(n2_o), n1_o + n2_o, df$n_o)
    # n_r_resolved <- ifelse(!is.na(n1_r) & !is.na(n2_r), n1_r + n2_r, df$n_r)
    #
    # if (score_mode) {
    #   p_o <- if ("orig_p_value"  %in% names(df)) df$orig_p_value  else rep(NA_real_, nrow(df))
    #   p_r <- if ("repli_p_value" %in% names(df)) df$repli_p_value else rep(NA_real_, nrow(df))
    #
    #   pScept <- mapply(compute_pScept_cor,
    #                    es_o = df$es_o, se_o = df$se_o, es_r = df$es_r, se_r = df$se_r,
    #                    n_o = n_o_resolved, n_r = n_r_resolved,
    #                    p_o = p_o, p_r = p_r)
    #
    #   skep_level <- ReplicationSuccess::levelSceptical(level = 0.05,
    #                                                    alternative = "two.sided",
    #                                                    type = "golden")
    #
    #   pattern_met <- if ("repli_pattern_criteria_met" %in% names(df)) {
    #     df$repli_pattern_criteria_met
    #   } else {
    #     sign(df$es_o) == sign(df$es_r)
    #   }
    #
    #   success <- !is.na(pScept) & !is.na(pattern_met) & (pScept < skep_level) & pattern_met
    #   na_mask <- is.na(pScept) | is.na(pattern_met)
    #
    # } else {
    #   pScept <- mapply(compute_pScept_z,
    #                    es_o = df$es_o, se_o = df$se_o, es_r = df$es_r, se_r = df$se_r,
    #                    n_o = n_o_resolved, n_r = n_r_resolved)
    #
    #   skep_level <- ReplicationSuccess::levelSceptical(level = alpha,
    #                                                    alternative = "one.sided",
    #                                                    type = "golden")
    #
    #   success <- !is.na(pScept) & (pScept <= skep_level)
    #   na_mask <- is.na(pScept)
    # }
    #
    # outcome_detailed <- ifelse(na_mask, NA_character_,
    #                            ifelse(success,
    #                                   "Evidence for replication success",
    #                                   "No evidence for replication success"))
    # outcome <- ifelse(na_mask, NA_character_,
    #                   ifelse(success, "Replication success", "Replication failure"))
    #
    # data.frame(outcome=outcome, outcome_detailed=outcome_detailed, quantitative_value = pScept)
    # }

  # Snapshot hybrid ----

  assess_snaphybrid <- function(df, threshold=0.75, pr_0 = 0.25, pr_s = 0.25, pr_m = 0.25){
    # Calculate prior for large effect based on other priors
    pr_l <- 1 - (pr_0 + pr_s + pr_m)
    stopifnot(pr_0 + pr_s + pr_m + pr_l == 1)

    effect_labels <- c("null effect", "small effect", "medium effect", "large effect")

    base_eligible <- if ("snaphybrid_base_eligible" %in% names(df)) {
      df$snaphybrid_base_eligible
    } else {
      rep(TRUE, nrow(df))
    }

    pSnap   <- rep(NA_real_, nrow(df))
    snap_idx <- rep(NA_integer_, nrow(df))

    for (i in seq_len(nrow(df))) {

      ## Match the 225-pair analysis sample used by Van Aert
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

      ## Snapshot uses Fisher-z correlations and raw sample-size variances
      es <- atanh(c(r_o, r_r))
      se <- 1 / sqrt(c(n_o, n_r) - 3)

      ## Orient the original positively and reverse the replication with it
      if (es[1] < 0) {
        es <- -es
      }

      ## Van Aert applies snapshot only when the original is significant:
      ## one-sided p < .025 after orientation
      p_original <- stats::pnorm(
        es[1] / se[1],
        lower.tail = FALSE
      )

      if (!is.finite(p_original) || p_original >= 0.025) {
        next
      }

      vals <- tryCatch(
        hybrid_snap(
          es = es,
          se = se,
          alpha = 0.025
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

  assess_Likratio <- function(df, lr = 8){

    scale_here <- unique(df$scale)
    to_r    <- function(es) if (identical(scale_here, "fisherz")) tanh(es) else es
    r_to_d  <- function(r)  2 * r / sqrt(1 - r^2)

    compute_Likratio <- function(es_o, es_r, n_o, n_r, lr) {
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
      # so the warning is suppressed.
      Dtie_res <- suppressWarnings(
        tryCatch(Dtie(n_o = n_o, lr = lr), error = function(e) NA_real_)
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
                     MoreArgs = list(lr = lr))

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
    # (repli_binary.R: "a shift in corr. interpretation"). The test must
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

    ### Success criterion, matching repli_binary.R's mapping exactly:
    ### Equivalence / Trivial Difference -> success; Difference -> failure;
    ### Indeterminacy -> NA (neither test conclusive)
    outcome <- ifelse(is.na(ct), NA_character_,
                      ifelse(ct %in% c("Equivalence", "Trivial Difference"), "Replication success",
                             ifelse(ct == "Difference", "Replication failure",
                                    ifelse(ct == "Indeterminacy", NA_character_, NA_character_))))

    data.frame(outcome = outcome, outcome_detailed = ct,
               quantitative_value = p_equiv)
  }

  # Minimum-effect test ----
  assess_minimum_effect <- function(df, h = NULL, shrinkage=0.25, alpha = 0.025){

  # Minimum-effect testing is the mirror image of equivalence testing (TOST):
  # Equivalently, H0 is rejected when the (1 - 2*alpha) CI for the difference
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

  results <- lapply(metric, function(m) {
    res <- metric_fns[[m]](df)
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

  bf_metrics <- c("BF_JZS", "BFe", "BFr", "BFs", "snaphybrid")
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

