library(mvnfast)
library(synlik)
library(latex2exp)
library(ggplot2)
library(cowplot)
library(nlme)

n <- 2000
truealpha <- 3
truebeta <- 1
truep <- 0.3

y.obs <- sample(c(0,1), n,
                prob=c(truep, 1-truep),replace=T) * rgamma(n,
                                                            truealpha,
                                                            truebeta)
S.obs <- c(mean(y.obs), sd(y.obs), mean(y.obs==0))

synthetic.likelihood <- function(pars, N, obs.length, S.obs) {
  Ss <- matrix(nrow=N, ncol=length(S.obs))
  
  for (i in 1:N) {
    xs <- sample(c(0,1),
                 obs.length,
                 prob=c(pars[3], 1-pars[3]), replace=T) * rgamma(n,
                                                                 pars[1],
                                                                 pars[2])
    Ss[i,] <- c(mean(xs), sd(xs), mean(xs==0))
  }
  
  mu <- colMeans(Ss)
  sigma <- cov(Ss)
  
  return(dmvn(S.obs, mu, sigma, log=TRUE))
}

likfun <- function(par) {
  return(synthetic.likelihood(par, 100, n, S.obs))
}

res <- synlik:::ml(likFun = likfun,
                   initPar = c(5, 5, 0.5), # Initial parameters
                   initCov = diag(c(1,1,0.1))^2, # Covariance matrix of the proposal
                   np = 50, # Number of simulated parameters
                   niter = 150) # Max num of iterations

df <- data.frame(res$estim)
colnames(df) <- c('alpha', 'beta', 'p')
df['Iteration'] = 1:nrow(df)

plot_1 <- ggplot(df, aes_string(x = 'Iteration', y = 'alpha')) +
  geom_point() +
  ylab(TeX('$\\alpha$')) +
  geom_hline(yintercept=truealpha, color = "red") +
  theme_half_open() +
  theme(text = element_text(size = 20), axis.title.y = element_text(angle = 0, size = 30))

plot_2 <- ggplot(df, aes_string(x = 'Iteration', y = 'beta')) +
  geom_point() +
  ylab(TeX('$\\beta$')) +
  geom_hline(yintercept=truebeta, color = "red") +
  theme_half_open() +
  theme(text = element_text(size = 20), axis.title.y = element_text(angle = 0, size = 30))

plot_3 <- ggplot(df, aes_string(x = 'Iteration', y = 'p')) +
  geom_point() +
  ylab(TeX('$p$')) +
  geom_hline(yintercept=truep, color = "red") +
  theme_half_open() +
  theme(text = element_text(size = 20), axis.title.y = element_text(angle = 0, size = 30))

plot <- plot_grid(plot_1, plot_2, plot_3, ncol = 3, nrow = 1)
tmpfile <- tempfile("no_robust_optimisation",
                    tmpdir = "/Users/jake/Documents/temp/dissertation/report/images/sl/gamma_example",
                    fileext = ".pdf")
save_plot(tmpfile, plot, base_height=4, base_width=12)

loss.fun <- function(x, s) {
  out <- sum(dnorm(x, s[1], s[2], log=TRUE)) - 0.5*sum((s[3] - as.numeric(x == 0))^2)
  return(-out)
}

likfun <- function(param) {
  par <- param
  
  sim <- sample(c(0,1), n, prob=c(par[3], 1-par[3]), replace=TRUE) * rgamma(n, par[1], par[2])
  s.sim <- c(mean(sim), sd(sim), mean(sim==0))
  
  H <- fdHess(s.sim, loss.fun, x = sim)$Hessian
  
  gr <- cbind(-sim/(s.sim[2]^2) + s.sim[1]/(s.sim[2]^2),
              1/s.sim[2] - (sim - s.sim[1])^2/(s.sim[2]^3),
              -as.numeric(sim==0) + sim[3])
  
  V <- cov(gr)*n
  SIG <- solve(H %*% solve(V, H))
  
  out <- dmvn(S.obs, s.sim, SIG, log = TRUE)
  
  return(out)
}

res <- synlik:::ml(likFun = likfun,
                   initPar = c(5, 5, 0.5), # Initial parameters
                   initCov = diag(c(1,1,0.1))^2, # Covariance matrix of the proposal
                   np = 50, # Number of simulated parameters
                   niter = 150) # Max num of iterations

df <- data.frame(res$estim)
colnames(df) <- c('alpha', 'beta', 'p')
df['Iteration'] = 1:nrow(df)

plot_1 <- ggplot(df, aes_string(x = 'Iteration', y = 'alpha')) +
  geom_point() +
  ylab(TeX('$\\alpha$')) +
  geom_hline(yintercept=truealpha, color = "red") +
  theme_half_open() +
  theme(text = element_text(size = 20), axis.title.y = element_text(angle = 0, size = 30))

plot_2 <- ggplot(df, aes_string(x = 'Iteration', y = 'beta')) +
  geom_point() +
  ylab(TeX('$\\beta$')) +
  geom_hline(yintercept=truebeta, color = "red") +
  theme_half_open() +
  theme(text = element_text(size = 20), axis.title.y = element_text(angle = 0, size = 30))

plot_3 <- ggplot(df, aes_string(x = 'Iteration', y = 'p')) +
  geom_point() +
  ylab(TeX('$p$')) +
  geom_hline(yintercept=truep, color = "red") +
  theme_half_open() +
  theme(text = element_text(size = 20), axis.title.y = element_text(angle = 0, size = 30))


plot <- plot_grid(plot_1, plot_2, plot_3, ncol = 3, nrow = 1)
tmpfile <- tempfile("robust_optimisation",
                    tmpdir = "/Users/jake/Documents/temp/dissertation/report/images/sl/gamma_example",
                    fileext = ".pdf")
save_plot(tmpfile, plot, base_height=4, base_width=12)