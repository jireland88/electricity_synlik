library(MASS)

# Function to return the trajectory values from xs that are on given weekdays
extract_days <- function(xs, day) {
  elements <- seq(day+1, length(xs), 7)
  return(xs[elements])
  
  #ys <- vector()
  #  for (i in seq_along(xs)) {
  #      if ((i - 1) %% 7 %in% days) {
  #          ys <- c(ys, xs[i])
  #      }
  #  }
  #  return(ys)
}

# Calculates the mean ratio between day t and day t+step
mean_rev <- function(xs, step) {
    ys <- vector()
    for (i in seq_len(length(xs) - step)) {
        ys <- c(ys, xs[i + step] / xs[i])
    }
    return(mean(ys))
}

polynomial_regression <- function(ys, xs, degree) {
  n <- length(xs)
  
  X <- matrix(0, n, degree)
  
  X <- cbind(rep(1, n), X)
  
  for (i in 1:degree) {
    X[, i+1] = xs^i
  }
  
  beta <- solve(t(X) %*% X) %*% t(X) %*%  ys
  
  return(beta)
}

loss_poly <- function(ys, xs, beta) {
  n <- length(xs)
  degree <- length(beta) - 1
  
  X <- matrix(0, n, degree)
  
  X <- cbind(rep(1, n), X)
  
  for (i in 1:degree) {
    X[, i+1] = xs^i
  }
  
  ys_hat <- X %*% beta
  
  loss <- sum((ys - ys_hat)^2)
  
  return(loss)
}

gradient_poly <- function(y, x, beta) {
  degree <- length(beta) - 1
  x_new <- c(1)
  
  for (i in 1:degree) {
    x_new <- c(x_new, x^i)
  }
  
  grad <- (-2 * t(x_new) * y) + (2 * t(x_new) * (x_new %*% beta)[1, 1])
  
  return(grad)
}

# Computes the linear regression:
# s(t) = beta_0 + beta_1 s(t-1) + ... + beta_shift s(t-shift) + error
autoregress <- function(obs, shift, force_intercept_origin = FALSE) {
    n <- length(obs)

    start = 1
    X <- matrix(0, n - shift, shift)

    if (force_intercept_origin == FALSE) {
      start = 2
      X <- cbind(rep(1, n - shift), X)    
    }

    for (i in start:(shift + 1)) {
        X[, i] <- obs[(i - 1):(n - shift + (i - 2))]
    }

    y <- obs[(shift + 1):n]

    beta <- solve(t(X) %*% X, t(X) %*%  y)

    return(beta)
}

# Computes the loss function for the autoregression above
loss_autoregress <- function(obs, shift, beta, force_intercept_origin = FALSE) {
    n <- length(obs)

    start = 1
    X <- matrix(0, n - shift, shift)

    if (force_intercept_origin == FALSE) {
      start = 2
      X <- cbind(rep(1, n - shift), X)    
    }

    for (i in start:(shift + 1)) {
        X[, i] <- obs[(i - 1):(n - shift + (i - 2))]
    }

    y <- obs[(shift + 1):n]
    y_hat <- X %*% beta

    loss <- sum((y - y_hat)^2)

    return(loss)
}

# Computes the gradient for the autoregression above
gradient_autoregress <- function(x, y, beta, force_intercept_origin = FALSE) {
    if (force_intercept_origin == FALSE) {
      x <- c(1, x)
    }

    grad <- (-2 * t(x) * y) + (2 * t(x) * (x %*% beta)[1, 1])

    return(grad)
}

# Get rid of values <= 0 and fit a gamma distribution on what's left
# getting rid of <= 0 values is now redundant I as fit the exponentiated trajectory
gamma <- function(xs, approx = FALSE) {
  #xs_gam <- xs[xs > 0]

  s <- log(mean(xs)) - mean(log(xs))
  k <- ( 3 - s + sqrt(((s - 3) ^ 2) + (24 * s)) ) / ( 12 * s )
  
  theta = mean(xs) / k
  
  par <- c(k, 1 / theta)
  
  if (!approx) {
    par <- fitdistr(xs, "gamma",
                    start = list(shape = k, rate = 1/theta))$estimate
  }
  
  return(par)
}

# Compute loss function for gamma above
loss_gamma <- function(xs, shape, rate) {
    xs_gam <- xs[xs > 0]

    return(-sum(dgamma(xs_gam, shape, rate, log = TRUE)))
}

# Compute gradient for gamma above
gradient_gamma <- function(x, shape, rate) {
    return(-c(log(rate * x) - digamma(shape), (shape / rate) - x))
}

# Compute loss function for a normal distribution
loss_normal <- function(xs, mu, sigma) {
  return(-sum(dnorm(xs, mu, sigma, log = TRUE)))
}

# Compute gradient for a normal distribution
gradient_normal <- function(x, mu, sigma) {
  grad <- c(-x / (sigma ^ 2) + mu / (sigma ^ 2),
            (1 / sigma) - ((x - mu) ^ 2) / (sigma ^ 3))

  return(grad)
}

# Compute the quadratic loss function
loss_quadratic <- function(xs, par) {
  return(sum((xs - par) ^ 2) / 2)
}

# Compute the gradient for the quadratic loss function
gradient_quadratic <- function(x, par) {
  return(par - x)
}

# return the ordered differences, that is y_{i+1} - y_i ordered from smallest to largest
ordered_differences <- function(xs) {
  diffs <- diff(xs)
  diffs_ordered <- diffs[order(diffs)]
  return(diffs_ordered)
}

spike_series <- function(xs) {
  d <- abs(diff(xs))
  spike_series <- c(0, d) + c(d, 0)
  return(spike_series)
}

# apologies for the name. This is the p_hat statistic in the report.
spike_test <- function(xs, n_std) {
  spike_series <- spike_series(xs)
  mu = mean(spike_series)
  sigma = sd(spike_series)
  n_spikes <- sum(spike_series > mu + n_std*sigma)
  
  return(n_spikes)
}