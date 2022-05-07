library(mvnfast)
library(latex2exp)
library(pheatmap)
library(RColorBrewer)
library(grid)
library(cowplot)

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

n <- 2000
truealpha <- 3
truebeta <- 1
truep <- 0.3

model_func <- function(pars) {
  pars[1] <- exp(pars[1])
  pars[2] <- exp(pars[2])
  pars[3] <- exp(pars[3]) / (1 + exp(pars[3]))

  obs <- sample(c(0, 1), n, c(pars[3], 1-pars[3]), replace = T) * rgamma(n, pars[1], pars[2])
  return(obs)
}

stats <- function(xs) {
  return(c(mean(xs), sd(xs), mean(xs==0)))
}

true_par <- c(log(truealpha), log(truebeta), log(truep/(1-truep)))

n_stats <- 3

cor_sp <- stat_param_cor(true_par, diag(c(1, 1, 1))^2,
                         stats, 10000, n_stats, model_func)

cor_ss <- stat_cor(true_par, diag(c(1, 1, 1))^2, stats,
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

cn <- c(' $\\alpha$ ', ' $\\beta$ ', ' $p$ ')
rn <- c('$\\bar{x}(\\cdot)$', 'sd($\\cdot$)', '$p$ est')

palette <- colorRampPalette(brewer.pal(7, name = 'Blues'))(100)

h_ss <- pheatmap(cor_ss, cluster_cols = FALSE, cluster_rows = FALSE,
                 main='Abs Correlation: statistics vs statistics',
                 fontsize=20, fontsize_row = 40, fontsize_col = 40,
                 color=palette, labels_row = TeX(rn), labels_col = TeX(rn))

h_sp <- pheatmap(cor_sp, cluster_cols = FALSE, cluster_rows = FALSE,
                 main='Abs Correlation: statistics vs parameters',
                 fontsize=20, fontsize_row = 40, fontsize_col = 40,
                 color=palette, labels_row = TeX(rn), labels_col = TeX(cn))

save_pheatmap_pdf <- function(x, filename, width=12, height=12) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_pdf(h_sp, "../report/images/sl/gamma_example/cor_sp.pdf")
save_pheatmap_pdf(h_ss, "../report/images/sl/gamma_example/cor_ss.pdf")

