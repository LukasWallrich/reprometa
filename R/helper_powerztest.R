# This function was copied from the source code of the R package "biostatUZH"
power.z.test <- function (n = NULL,
                          delta = NULL,
                          sd = 1,
                          sig.level = 0.05,
                          power = NULL,
                          type = c("two.sample", "one.sample"),
                          alternative = c("two.sided", "one.sided")
){

  # Check exactly one of n, delta, sd, power, sig.level is NULL
  # Adapted from stats::power.t.test
  if (sum(vapply(list(n, delta, sd, power, sig.level), is.null,
                 logical(1L))) != 1)
    stop("exactly one of 'n', 'delta', 'sd', 'power', and 'sig.level' must be NULL")
  # Check that n is larger than 1
  if(!is.null(n) && n < 1)
    stop("n must be larger than 1")
  # Check that sd is positive
  if(!is.null(sd) && sd <= 0)
    stop("sd must be positive")

  # Check that sig.level and power are either NULL or between 0 and 1
  assert_NULL_or_prob(sig.level)
  assert_NULL_or_prob(power)

  # set alternative and type for now
  type <- match.arg(type)
  tsample <- switch(type, one.sample = 1, two.sample = 2)# , paired = 1)
  alternative <- match.arg(alternative)
  tside <- switch(alternative, one.sided = 1, two.sided = 2)

  if(is.null(n)){
    n <- tsample * ((qnorm(power) + qnorm(1 - sig.level/tside)) * sd / delta)^2
  } else if(is.null(delta)){
    delta <- sd * sqrt(tsample / n) * (qnorm(power) + qnorm(1 - sig.level/tside))
  } else if(is.null(sd)){
    sd <- abs(delta) * sqrt(n) / (sqrt(tsample) * (qnorm(power) + qnorm(1 - sig.level/tside)))
  } else if(is.null(sig.level)){
    sig.level <- tside * (1 - pnorm(sqrt(n/tsample) * abs(delta)/sd - qnorm(power)))
  } else if(is.null(power)){
    power <- pnorm(sqrt(n/tsample) * abs(delta) / sd - qnorm(1 - sig.level/tside))
  }

  # Create notes (copied from stats::power.t.test)
  NOTE <- switch(type,
                 paired = "n is number of *pairs*, sd is std.dev. of *differences* within pairs",
                 two.sample = "n is number in *each* group",
                 NULL)
  METHOD <- paste(switch(type,
                         one.sample = "One-sample",
                         two.sample = "Two-sample",
                         paired = "Paired"),
                  "z test power calculation")

  # Return object (copied from stats::power.t.test)
  structure(list(n = n, delta = delta, sd = sd, sig.level = sig.level,
                 power = power, alternative = alternative, note = NOTE,
                 method = METHOD), class = "power.htest")

}

assert_NULL_or_prob <- function(x){
  if(!is.null(x) && (!is.numeric(x) || any(0 > x | x > 1)))
    stop(gettextf("'%s' must be numeric in [0, 1]", deparse1(substitute(x))),
         call. = FALSE)
  invisible(x)
}
