library(mvnfast)
library(latex2exp)
library(pheatmap)
library(RColorBrewer)
library(grid)
library(cowplot)
library(moments)
library(MASS)

source('models.r')
source('helper_functions.r')

stat_param_cor <- function(init_par, covar, stats, n_sim, n_stats, model_func) {
  pars <- rmvn(n_sim, init_par, covar)
  stats_vecs <- matrix(NA, n_sim, n_stats)
  cors <- matrix(NA, n_stats, length(init_par))
  
  for (i in seq(n_sim)) {
    stats_vecs[i, ] <- stats(model_func(pars[i, ]))
  }
  
  pars <- pars[!is.na(rowSums(stats_vecs)),]
  stats_vecs <- stats_vecs[!is.na(rowSums(stats_vecs)),]
  
  for (p in seq_along(init_par)) {
    for (s in seq(n_stats)) {
      cors[s, p] = abs(cor(pars[, p], stats_vecs[, s], use = "complete.obs"))
    }
  }
  
  return(cors)
}

stat_cor <- function(init_par, covar, stats, n_sim, n_stats, model_func) {
  pars <- rmvn(n_sim, init_par, covar)
  stats_vecs <- matrix(NA, n_sim, n_stats)
  
  for (i in seq(n_sim)) {
    stats_vecs[i, ] <- stats(model_func(pars[i, ]))
  }
  
  stats_vecs <- stats_vecs[!is.na(rowSums(stats_vecs)),]
  
  c <- abs(cor(stats_vecs, use = "complete.obs"))
  
  return(c)
}

obs_length <- 1000

true_par <- c(2.852, -0.089, -0.192, log(0.112 / (1 - 0.112)), log(0.114),
              0.103, log(0.542), log(0.313 / (1 - 0.313)), log(0.453),
              log(0.95 / (1 - 0.95)))

obs <- do.call(regime_jump_model, as.list(append(true_par, 1000)))
obs_ord <- cbind(rep(1, length(obs)-1),
                 poly(ordered_differences(obs), 3))

model_func <- function(pars) {
  obs <- do.call(regime_jump_model, as.list(append(pars, 1000)))
  return(obs)
}

stats <- function(xs) {
  sat <- extract_days(xs, c(5))
  sun <- extract_days(xs, c(6))
  
  mu <- mean(xs)
  
  diffs_1 <- diff(xs)
  diffs_2 <- diff(diffs_1)
  
  X <- cbind(rep(1, obs_length - 3),
             xs[1:(obs_length-3)],
             xs[2:(obs_length-2)],
             xs[3:(obs_length-1)])
  
  reg <- lm.fit(X, xs[4:obs_length])$coefficients
  
  shape <- lm.fit(obs_ord, ordered_differences(xs))$coefficients
  
  return(c(mu, mean(sat) - mu, mean(sun) - mu, sd(xs), IQR(xs),
           max(abs(xs)), mean(diffs_1), max(abs(diffs_1)), sd(diffs_1),
           skewness(diffs_1), kurtosis(diffs_1), mean(diffs_2),
           max(abs(diffs_2)), sd(diffs_2), skewness(diffs_2),
           kurtosis(diffs_2), reg, as.vector(acf(xs, 3)$acf),
           skewness(xs), mean(xs) - median(xs), gamma(exp(xs), TRUE),
           spike_test(xs, 2), shape))
}

n_stats = 33

cor_sp <- stat_param_cor(true_par, diag(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1))^2,
                         stats, 1000, n_stats, model_func)

cor_ss <- stat_cor(true_par, diag(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1))^2,
                   stats, 1000, n_stats, model_func)

draw_colnames_0 <- function (coln, gaps, ...) {
  coord <- pheatmap:::find_coordinates(length(coln), gaps)
  x     <- coord$coord - 0.5 * coord$size
  res   <- grid::textGrob(
    coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"),
    vjust = 1.25, rot = 0, gp = grid::gpar(...)
  )
  return(res)
}

cn <- c('$\\mu_0$', '$\\beta_1$', '$\\beta_2$', '$\\alpha_0$', '$\\sigma_0$',
        '$\\mu_1$', '$\\sigma_1$', '$\\alpha_{-1}$', '$\\sigma_{-1}$', '$p$')

rn <- c('$\\bar{x}(\\mathbf{y})$', '$\\bar{x}_5(\\mathbf{y}) - \\bar{x}(\\mathbf{y})$',
        '$\\bar{x}_6(\\mathbf{y}) - \\bar{x}(\\mathbf{y})$', 'sd($\\mathbf{y}$)',
        'IQR($\\mathbf{y}$)', 'max($|\\mathbf{y}|$)',
        '$\\bar{x}(\\mathbf{d}^1$)', 'max($|\\mathbf{d}^1|$)',
        'sd($\\mathbf{d}^1$)', '$sk($\\mathbf{d}^1$)',
        'k($\\mathbf{d}^1$)', '$\\bar{x}(\\mathbf{d}^2)$',
        'max($|\\mathbf{d}^1|$)', 'sd($\\mathbf{d}^2$)',
        '$sk($\\mathbf{d}^2$)', 'k($\\mathbf{d}^2$)',
        '$\\xi_0(\\mathbf{y})$', '$\\xi_1(\\mathbf{y})$', '$\\xi_2(\\mathbf{y})$',
        '$\\xi_3(\\mathbf{y})$', '$\\nu_0(\\mathbf{y})$', '$\\nu_1(\\mathbf{y})$',
        '$\\nu_2(\\mathbf{y})$', '$\\nu_3(\\mathbf{y})$', 'sk($\\mathbf{y}$)',
        '$\\bar{x}(\\mathbf{y}) - m(\\mathbf{y})$',
        '$\\theta_1(\\mathbf{y})$', '$\\mathbf{\\theta}_2(\\mathbf{y})$',
        '$\\hat{p}(\\mathbf{y})$', '$\\gamma_0(\\mathbf{y})$', '$\\gamma_1(\\mathbf{y})$',
        '$\\gamma_2(\\mathbf{y})$', '$\\gamma_3(\\mathbf{y})$')

palette <- colorRampPalette(brewer.pal(7, name = 'Blues'))(100)

.draw_col <- get("draw_colnames", asNamespace('pheatmap'))

assignInNamespace(
  x = "draw_colnames",
  value = "draw_colnames_0",
  ns = asNamespace("pheatmap")
)

h_sp <- pheatmap(cor_sp, cluster_cols = FALSE, cluster_rows = FALSE,
                 main='Abs Correlation: statistics vs parameters',
                 labels_col = TeX(cn), labels_row = TeX(rn),
                 fontsize=20, fontsize_row = 20, fontsize_col = 20,
                 color=palette)

assignInNamespace(
  x = "draw_colnames",
  value = ".draw_col",
  ns = asNamespace("pheatmap")
)

h_ss <- pheatmap(cor_ss, cluster_cols = FALSE, cluster_rows = FALSE,
                 main='Abs Correlation: statistics vs statistics',
                 fontsize=20, fontsize_row = 20, fontsize_col = 20,
                 color=palette, labels_row = TeX(rn), labels_col = TeX(rn))

pdf("../report/images/fitting/full_model/cor_ss.pdf", width=12, height=12)
h_ss
dev.off()

pdf("../report/images/fitting/full_model/cor_sp.pdf", width=12, height=12)
h_sp
dev.off()


