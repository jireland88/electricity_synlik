library(reshape2)
library(ggplot2)
library(cowplot)

ss_full <- read.table('ss_full.txt')

ss_full <- sweep(ss_full, 2, true_par)

true_par <- c(2.852, -0.089, -0.192, log(0.112 / (1 - 0.112)), log(0.114),
              log(0.103), log(0.542), log(0.313 / (1 - 0.313)), log(0.453),
              log(0.95 / (1 - 0.95)))

options(repr.plot.width=10, repr.plot.height=10)
p <- ggplot(melt(ss_full), aes(x=variable, y=value)) +
  geom_boxplot() +
  coord_flip() + 
  geom_hline(yintercept=0, col='Red') +
  theme_half_open() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "Full Model Simulation Study", x = "", y = "") +
  scale_x_discrete(labels=c(expression(mu[0]),
                            expression(beta[1]),
                            expression(beta[2]),
                            expression(alpha[0]),
                            expression(sigma[0]),
                            expression(mu[1]),
                            expression(sigma[1]),
                            expression(alpha[-1]),
                            expression(sigma[-1]),
                            expression(p)))

pdf("../report/images/fitting/full_model/ss.pdf", width=12, height=12)
p
dev.off()
