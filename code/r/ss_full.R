library(mvnfast)
library(nlme)
library(synlik)
library(moments)

source("models.r")
source("helper_functions.r")

obs_length <- 100
true_par <- c(2.852, -0.089, -0.192, log(0.112 / (1 - 0.112)), log(0.114),
              log(0.103), log(0.542), log(0.313 / (1 - 0.313)), log(0.453),
              log(0.95 / (1 - 0.95)))

num_sims <- 50

ss <- matrix(nrow=num_sims, ncol=length(true_par))

for (j in 1:num_sims) {
  print(j)
  obs <- do.call(regime_jump_model, as.list(append(true_par, obs_length)))
  
  obs_ord_poly <- cbind(rep(1, length(obs)-1),
                        poly(ordered_differences(obs), 3))
  
  stats <- function(xs) {
    sat <- extract_days(xs, 5)
    sun <- extract_days(xs, 6)
    
    mu <- mean(xs)
    sigma <- sd(xs)
    
    med <- median(xs)
    
    # autoregression
    
    X <- cbind(rep(1, obs_length - 3),
               xs[1:(obs_length-3)],
               xs[2:(obs_length-2)],
               xs[3:(obs_length-1)])
    
    reg <- lm.fit(X, xs[4:obs_length])$coefficients
    
    gam <- gamma(exp(xs), TRUE)
    
    #shape <- lm.fit(obs_od_poly, ordered_differences(xs))$coefficients
    
    shape <- lm.fit(obs_ord_poly, ordered_differences(xs))$coefficients
    
    #spike_series_xs <- spike_series(xs)
    #spike_reg <- lm.fit(spike_series_obs_ord,
    #                    spike_series_xs[order(spike_series_xs)])$coefficients
    
    #dif_reg <- lm.fit(obs_ord_diff, ordered_differences(xs))$coefficients
    #spike <- length(spike_detection(xs, 0.45))
    
    ttr <- times_to_revert(xs)
    
    diffs <- diff(xs)
    
    # dropping: kurt, sd
    
    return(c(mu, mean(sat) - mu, mean(sun) - mu, reg[2],
             skewness(xs), reg[4], IQR(xs), max(abs(diffs)),
             max(abs(xs)), spike_test(xs, 2), shape))
  }
  
  synthetic_likelihood <- function(pars, N, obs_length, s_obs) {
    stats_matrix <- matrix(0, nrow=N, ncol=length(s_obs))
    
    for (i in 1:N) {
      sim <- do.call(regime_jump_model, as.list(append(pars, obs_length)))
      stats_matrix[i, ] = stats(sim)
    }
    
    mu <- colMeans(stats_matrix)
    Sigma <- cov(stats_matrix)
    
    out <- dmvn(s_obs, mu, Sigma, log = TRUE)
    
    return(out)
  }
  
  s_obs <- stats(obs)
  print(s_obs)
  
  N=100
  
  likfun <- function(params) {
    return(synthetic_likelihood(params, N, obs_length, s_obs))
  }
  
  print(likfun(true_par))
  
  res <- synlik:::ml(likFun = likfun,
                     initPar = true_par,
                     initCov = diag(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1))^2,
                     np = 50,
                     niter = 100)
  
  ss[j, ] <- res$estim[100, ]
}

write.table(ss, 'ss_full.txt')



