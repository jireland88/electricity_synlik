
library(latex2exp)
library(ggplot2)
library(cowplot)

regime.jump.model.constrained <- function(mu.0, beta.1, beta.2, alpha.0, sigma.0,
                                          mu.1, sigma.1, alpha.minus1, sigma.minus1,
                                          p, t.max) {
  
  ys <- rep(NA, t.max)
  
  regime <- 0
  x <- 0
  
  for (t in 1:t.max) {
    
    # deterministic part
    f <- mu.0
    day <- t %% 7
    
    if (day == 5) {
      # Saturday
      f <- f + beta.1
    } else if (day == 6) {
      # Sunday
      f <- f + beta.2
    }
    
    # stochastic part (select regime based on last regime, then calculate dx)
    dx <- 0
    if (regime == 0) {
      regime <- sample(c(0,1), size=1, replace=TRUE, prob=c(p, 1-p))
    } else if (regime == 1) {
      regime <- -1
    } else if (regime == -1) {
      regime <- 0
    }
    
    if (regime == 0) {
      dx <- (-alpha.0 * x) + (sigma.0 * rnorm(1))
    } else if (regime == 1) {
      dx <- mu.1 + (sigma.1 * rnorm(1))
    } else if (regime == -1) {
      dx <- (-alpha.minus1 * x) + (sigma.minus1 * rnorm(1))
    }
    
    x <- x + dx
    
    ys[t] <- f + x
  }
  
  return(ys)
}

xs <- 1:100
ys <- regime.jump.model.constrained(2.852, -0.089, -0.192, 0, 0.114, 0, 0, 0, 0, 1, 100)

df <- data.frame(x = xs, y = ys)

plot_1 <- ggplot(data=df, aes(x=x, y=y)) +
  geom_line() +
  ylab(TeX('$s(t)$')) +
  xlab(TeX('$t$')) +
  ggtitle(TeX('$p = 1, \\alpha_0 = 0$')) +
  theme_half_open() +
  theme(text = element_text(size = 20), axis.title.y = element_text(angle = 0, size = 20))

xs <- 1:100
ys <- regime.jump.model.constrained(2.852, -0.089, -0.192, 0.112, 0.114, 0, 0, 0, 0, 1, 100)

df <- data.frame(x = xs, y = ys)

plot_2 <- ggplot(data=df, aes(x=x, y=y)) +
  geom_line() +
  ylab(TeX('$s(t)$')) +
  xlab(TeX('$t$')) +
  ggtitle(TeX('$p = 1$')) +
  theme_half_open() +
  theme(text = element_text(size = 20), axis.title.y = element_text(angle = 0, size = 20))

plot <- plot_grid(plot_1, plot_2, ncol = 2, nrow = 1)
tmpfile <- tempfile("mean_reversion_example",
                    tmpdir = "/Users/jake/Documents/temp/dissertation/report/images/model",
                    fileext = ".pdf")
save_plot(tmpfile, plot, base_height=4, base_width=12)


xs <- 1:100
ys <- regime.jump.model.constrained(2.852, -0.089, -0.192, 0.112, 0.114, 0.103, 0.542, 0.313, 0.453, 0.95, 100)

df <- data.frame(x = xs, y = ys)

plot_1 <- ggplot(data=df, aes(x=x, y=y)) +
  geom_line() +
  ylab(TeX('$s(t)$')) +
  xlab(TeX('$t$')) +
  ggtitle('UK Parameters') +
  theme_half_open() +
  theme(text = element_text(size = 20), axis.title.y = element_text(angle = 0, size = 20))

ys <- regime.jump.model.constrained(2.852, -0.089, -0.192, 0.112, 0.114, 2, 0.542, 0.313, 0.453, 0.95, 100)

df <- data.frame(x = xs, y = ys)

plot_2 <- ggplot(data=df, aes(x=x, y=y)) +
  geom_line() +
  ylab(TeX('$s(t)$')) +
  xlab(TeX('$t$')) +
  ggtitle(TeX('$\\mu_0 = 2$')) +
  theme_half_open() +
  theme(text = element_text(size = 20), axis.title.y = element_text(angle = 0, size = 20))

ys <- regime.jump.model.constrained(2.852, -0.089, -0.192, 0.112, 0.114, 2, 0.542, 0.313, 0.453, 0.99, 100)

df <- data.frame(x = xs, y = ys)

plot_3 <- ggplot(data=df, aes(x=x, y=y)) +
  geom_line() +
  ylab(TeX('$s(t)$')) +
  xlab(TeX('$t$')) +
  ggtitle(TeX('$\\mu_0 = 2, p = 0.99$')) +
  theme_half_open() +
  theme(text = element_text(size = 20), axis.title.y = element_text(angle = 0, size = 20))

plot <- plot_grid(plot_1, plot_2, plot_3, ncol = 3, nrow = 1)
tmpfile <- tempfile("spike_example",
                    tmpdir = "/Users/jake/Documents/temp/dissertation/report/images/model",
                    fileext = ".pdf")
save_plot(tmpfile, plot, base_height=4, base_width=12)