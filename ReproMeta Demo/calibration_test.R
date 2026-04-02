# Calibration of various metrics to the overal T1E alpha squared
library(ReplicationSuccess)

data("RProjects")
test_df <- RProjects[RProjects$`project` == "Psychology", ]

# FE Meta-analysis ----
eps <- 0.000001
alpha <- 0.025

# Take real data as an example
es_o <- test_df$fiso[1]
es_r <- test_df$fisr[1]
se_o <- test_df$se_fiso[1]
se_r <- test_df$se_fisr[1]

# Stouffer's method ----
p1 <- 1 - seq(eps, 1-eps, by = eps)
z1 <- p2z(p1, alternative = "one.sided")
z2 <- sqrt(2) * qnorm(1-alpha^2) - z1
p2 <- z2p(z2, alternative = "one.sided") 

plot(p1,p2, type = "l", xlim = c(0,0.25), ylim = c(0,0.25))
points(x = z2p(es_o / se_o, alternative = "one.sided"), 
       y = z2p(es_r / se_r, alternative = "one.sided"))

1-pnorm((es_o / se_o + es_r / se_r) / sqrt(2))


set.seed(2026)
# simulation: ----

alpha <- 0.025
nsim <- 10^6
# Under the null hypothesis, both studies have a zero effect

theta_o <- rnorm(nsim, mean = 0, sd = 1)
theta_r <- rnorm(nsim, mean = 0, sd = 1)

# Calibrate the success rate to alpha^2
zcrit <- qnorm(1-alpha^2)

## Prediction interval ----
# Check the success rate that we would have at random: Prediction interval
success_PI <- zcrit <= (theta_r - theta_o) / sqrt(1^2 + 1^2)
mean(success_PI)

# Check the success rate that we would have at random: Original effect in Replication Confidence interval

success_CIr <- zcrit <= (theta_o - theta_r) / sqrt(1^2)
mean(success_CIr)

# Check the success rate that we would have at random: Replication effect in Original Confidence interval

success_CIo <- zcrit <= (theta_r - theta_o) / sqrt(1^2)
mean(success_CIo)

## Small telescopes ----
# Null hypothesis: is the replication study significantly different than 
# the effect that the original study would have had 33% power to detect
# H_0: theta_r <= theta_33
# H_1: theta_r > theta_33
set.seed(10)
sapply(c(0, 0.0001, 0.005, 0.0075, 0.01, 0.05, 0.1, 1), function(eps) {
# Under the null hypothesis, both studies have a zero effect
theta_33 <- eps*(qnorm(alpha, lower.tail = FALSE) + qnorm(1/3, lower.tail = FALSE)) # eps is equal to sd/sqrt(n_o) -> n_o is the 


success_SmallT <- qnorm(1-(alpha^2), mean = 0, sd = 1) <= (theta_r - theta_33) / 1
return(mean(success_SmallT))
})

library(parallel)

nsim <- 1000000
nrep <- 1000
n_cores <- detectCores() - 1

smalltelescopes_uniroot <- function(eps, theta_r, alpha = 0.025) {
  theta_33       <- eps * (qnorm(alpha, lower.tail = FALSE) + qnorm(1/3, lower.tail = FALSE))
  success_SmallT <- qnorm(1 - alpha^2) <= (theta_r - theta_33)
  return(mean(success_SmallT) - alpha^2)
}

store_mean <- mclapply(1:nrep, function(i) {
  set.seed(i)
  theta_r <- rnorm(nsim)
  uniroot(smalltelescopes_uniroot, interval = c(-1, 10), 
          tol = 1e-9, theta_r = theta_r)$root
}, mc.cores = n_cores)

mean(unlist(store_mean))
sd(unlist(store_mean))



## For small telescopes the miscalibration depends on theta_33. 
## When it's equal to zero it converges to alpha^2. When theta_33 goes to 1, overall T1E is zero 

## BF_JZS ----

set.seed(10)
theta_o <- rnorm(nsim, mean = 0, sd = 1)
theta_r <- rnorm(nsim, mean = 0, sd = 1)

# Computes integrand of the likelihood and the Cauchy prior (JZS prior)
integrate_jzs <- function(theta, es, se){
  stats::dnorm(x = es, mean = theta, sd = se) *
    stats::dcauchy(theta, location = 0, scale = 1)
}

# Divide by marginal likelihood under the null to obtain BF
compute_bfjzs <- function(es, se) {
  stats::integrate(integrate_jzs, lower=-Inf, upper=Inf, es=es, se=se)$value /
    stats::dnorm(x=es, mean=0, sd=se)
}

BFJZS_ori <- mapply(compute_bfjzs, es = theta_o, se = 1)
BFJZS_rep <- mapply(compute_bfjzs, es = theta_r, se = 1)
sapply(c(1, 2, 2.5, 3, 4, 5, 10), function(gamma) {
  mean(BFJZS_ori > gamma & BFJZS_rep > gamma)
})
# miscalibration factor for JZS: when gamma is set to 3, then overall T1E = 0.000266
# When using the standard threshold for replication success

JZS_uniroot <- function(gamma, alpha = 0.025) {
  mean(BFJZS_ori > gamma & BFJZS_rep > gamma) - alpha^2
}

result_BFJZS <- uniroot(JZS_uniroot, interval = c(1, 100), tol = 1e-6)
result_BFJZS$root
# Optimal gamma for overall T1E control at alpha^2 is 2.351871

## BFe ----

# For this metric we need an heterogeneity factor --> which one should we pick?
BFe <- BayesRep::BFe(to=theta_o, so=1, tr=theta_r, sr=1, tau=0.3)

sapply(c(1, 2, 2.5, 3, 5, 10, 20, 100), function(gamma){
  mean(BFe > gamma)
})
# For tau = 0.3, the miscalibrated overall T1E is 0

# Which threshold gamma fixes the overall T1E to alpha^2?
BFe_uniroot <- function(gamma, alpha = 0.025) {
  mean(BFe > gamma) - alpha^2
}

result_BFe <- uniroot(BFe_uniroot, interval = c(1, 100), tol = 1e-6)
result_BFe$root

BFe_1 <- BayesRep::BFe(to=theta_o, so=1, tr=theta_r, sr=1, tau=1)

sapply(c(1, 2, 2.5, 3, 5, 10, 20, 100), function(gamma){
  mean(BFe_1 > gamma)
})

# For tau = 1, the miscalibrated overall T1E is 1

# Which threshold gamma fixes the overall T1E to alpha^2?
BFe_1_uniroot <- function(gamma, alpha = 0.025) {
  mean(BFe_1 > gamma) - alpha^2
}

result_BFe_1 <- uniroot(BFe_1_uniroot, interval = c(0.1, 100), tol = 1e-9)
result_BFe_1$root

## BFr ----


BFr <- 1 / BayesRep::BFr(to=theta_o, so=1, tr=theta_r, sr=1)

sapply(c(1, 2, 2.5, 3, 5, 10, 20), function(gamma){
  mean(BFr > gamma)
})
# Gamma = 3 miscalibrates overall T1E to 0.029496


# Which threshold gamma fixes the overall T1E to alpha^2?
BFr_uniroot <- function(gamma, alpha = 0.025) {
  mean(BFr > gamma) - alpha^2
}

result_BFr <- uniroot(BFr_uniroot, interval = c(1, 100), tol = 1e-6)
result_BFr$root
# Ideal gamma 35.43926 (however this test does not consider two null hypotheses, and is not suggested to be used when original effect is zero)

## BFs ----

BFs <- 1 / BayesRep::BFs(to=theta_o, so=1, tr=theta_r, sr=1)
# Drop spectical BFs that did not converge (NaNs)
new_BFs <- BFs[!is.na(BFs)]

# Miscalibrated overall T1E value depends on gamma 
# (the threshold for success in the Bayesian framework)
sapply(c(1, 2, 2.5, 3, 5, 10, 20), function(gamma){
  mean(new_BFs > gamma)
})
# Miscalibrated BFs overall T1E when gamma=3 is : 0.0002286


# Which threshold gamma fixed the overall T1E to alpha^2
BFs_uniroot <- function(gamma, alpha = 0.025) {
  mean(new_BFs > gamma) - alpha^2
}

result_BFs <- uniroot(BFs_uniroot, interval = c(1, 100), tol = 1e-6)
result_BFs$root
# Calibrated BFs is 4.08
sd(new_BFs)
## snaphybrid ----
pr_0 <- pr_s <- pr_m <- 0.25
pr_l <- 1 - (pr_0 + pr_s + pr_m)

effect_labels <- c("null effect", "small effect", "medium effect", "large effect")

snap_results <- mapply(function(es_o, es_r, se_o, se_r) {
  vals  <- hybrid_snap(es = c(es_o, es_r), se = c(se_o, se_r))
  probs <- unlist(vals)
  idx   <- which.max(probs)
  c(pSnap = unname(probs[idx]), snap_idx = unname(idx))
}, theta_o, theta_r, 1, 1)

# snap_results: rows = pSnap/snap_idx, cols = studies
post_Snap <- snap_results[1, ]
# BF: posterior odds/ prior odds.
# prior odds: 0.25 / (1-0.25) = 1/3. -> Holds when flat prior on models
BF_Snap <- (post_Snap / (1 - post_Snap)) * ((1 - 0.25) / 0.25)

# Miscalibrated overall T1E value depends on gamma 
# (the threshold for success in the Bayesian framework)
sapply(c(1, 2, 2.5, 3, 5, 10, 20), function(gamma){
  mean(BF_Snap > gamma)
})
# Gamma = 3 miscalibrates the overall T1E rate to 0.097219
# Which threshold gamma fixed the overall T1E to alpha^2
BF_Snap_uniroot <- function(gamma, alpha = 0.025) {
  mean(BF_Snap > gamma) - alpha^2
}

result_BF_Snap <- uniroot(BF_Snap_uniroot, interval = c(1, 100), tol = 1e-6)
# BF For flat and equal priors: 4.7622 ideal gamma controlling for overall T1E
result_BF_Snap$root

## equivalence (?) ----
