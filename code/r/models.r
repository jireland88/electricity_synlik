regime_jump_model <- function(mu_0, beta_1, beta_2, alpha_0, sigma_0,
                              mu_1, sigma_1, alpha_neg_1, sigma_neg_1,
                              p, t_max) {

  # transform all constrained parameters
  mu_1 <- exp(mu_1)
  sigma_0 <- exp(sigma_0)
  sigma_1 <- exp(sigma_1)
  sigma_neg_1 <- exp(sigma_neg_1)
  alpha_0 <- exp(alpha_0) / (1 + exp(alpha_0))
  alpha_neg_1 <- exp(alpha_neg_1) / (1 + exp(alpha_neg_1))
  p <- exp(p) / (1 + exp(p))
  ys <- rep(NA, t_max)

  epsilons <- rnorm(t_max)

  # pre-calculate samples for efficiency
  enter_spike_samples <- sample(c(0, 1),
                                size = t_max,
                                replace = TRUE,
                                prob = c(p, 1 - p))
  
  regime <- 0
  x <- 0

  for (t in 1:t_max) {
    # deterministic part
    f <- mu_0
    day <- (t - 1) %% 7

    if (day == 5) {
      # Saturday
      f <- f + beta_1
    } else if (day == 6) {
      # Sunday
      f <- f + beta_2
    }

    # stochastic part (select regime based on last regime, then calculate dx)
    dx <- 0
    if (regime == 0) {
      regime <- enter_spike_samples[t]
    } else if (regime == 1) {
      regime <- -1
    } else if (regime == -1) {
      regime <- 0
    }

    if (regime == 0) {
      dx <- (-alpha_0 * x) + (sigma_0 * epsilons[t])
    } else if (regime == 1) {
      dx <- mu_1 + (sigma_1 * epsilons[t])
    } else if (regime == -1) {
      dx <- (-alpha_neg_1 * x) + (sigma_neg_1 * epsilons[t])
    }

    x <- x + dx

    ys[t] <- f + x
  }

  return(ys)
}

rjm_simple <- function(mu_0, beta_1, beta_2, alpha_0, sigma_0, t_max) {
  
  # transform all constrained parameters
  sigma_0 <- exp(sigma_0)
  alpha_0 <- exp(alpha_0) / (1 + exp(alpha_0))

  ys <- rep(NA, t_max)
  
  epsilons <- rnorm(t_max)
  
  x <- 0
  
  for (t in 1:t_max) {
    # deterministic part
    f <- mu_0
    day <- (t - 1) %% 7
    
    if (day == 5) {
      # Saturday
      f <- f + beta_1
    } else if (day == 6) {
      # Sunday
      f <- f + beta_2
    }
    
    # stochastic part
    dx <- (-alpha_0 * x) + (sigma_0 * epsilons[t])
    
    x <- x + dx
    
    ys[t] <- f + x
  }
  
  return(ys)
}