library(mvnfast)
library(latex2exp)
library(pheatmap)
library(RColorBrewer)
library(grid)
library(cowplot)

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

model_func <- function(pars) {
  obs <- do.call(rjm_simple, as.list(append(pars, 1000)))
  return(obs)
}

stats <- function(xs) {
  sat <- extract_days(xs, c(5))
  sun <- extract_days(xs, c(6))
  
  mu <- mean(xs)
  sigma <- sd(xs)
  
  
  reg <- autoregress(xs, 1)
  
  
  return(c(mu, mean(sat), mean(sun), reg, sigma))
}

true_par <- c(2.852, -0.089, -0.192, log(0.112 / (1 - 0.112)), log(0.114))

n_stats = 6

cor_sp <- stat_param_cor(true_par, diag(c(1, 1, 1, 1, 1))^2,
                         stats, 10000, n_stats, model_func)

cor_ss <- stat_cor(true_par, diag(c(1, 1, 1, 1, 1))^2, stats,
                   1000, n_stats, model_func)

draw_colnames_0 <- function (coln, gaps, ...) {
  coord <- pheatmap:::find_coordinates(length(coln), gaps)
  x     <- coord$coord - 0.5 * coord$size
  res   <- grid::textGrob(
    coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"),
    vjust = 1.25, rot = 0, gp = grid::gpar(...)
  )
  return(res)
}

assignInNamespace(
  x = "draw_colnames",
  value = "draw_colnames_0",
  ns = asNamespace("pheatmap")
)

cn <- c('$\\mu_0$', '$\\beta_1$', '$\\beta_2$', '$\\alpha_0$', '$\\sigma_0$')
rn <- c('$\\bar{x}(\\mathbf{y})$', '$\\bar{x}_5(\\mathbf{y})$', '$\\bar{x}_6(\\mathbf{y})$',
        '$\\xi_0(\\mathbf{y})$', '$\\xi_1(\\mathbf{y})$', 'sd($\\mathbf{y}$)')

palette <- colorRampPalette(brewer.pal(7, name = 'Blues'))(100)

h_ss <- pheatmap(cor_ss, cluster_cols = FALSE, cluster_rows = FALSE,
                 main='Abs Correlation: statistics vs statistics',
                 fontsize=30, fontsize_row = 40, fontsize_col = 40,
                 color=palette, labels_row = TeX(rn), labels_col = TeX(rn))

h_sp <- pheatmap(cor_sp, cluster_cols = FALSE, cluster_rows = FALSE,
                 main='Abs Correlation: statistics vs parameters',
                 labels_col = TeX(cn), labels_row = TeX(rn),
                 fontsize=30, fontsize_row = 40, fontsize_col = 40,
                 color=palette)

pdf("../report/images/fitting/simple_model/cor_ss.pdf", width=15, height=12)
h_ss
dev.off()

pdf("../report/images/fitting/simple_model/cor_sp.pdf", width=15, height=12)
h_sp
dev.off()


