library(GGally)
library(tidyverse)
library(here)
require(latex2exp)
library(KScorrect)
post_tibble <- read_csv('Posterior/log_adaptive_4.csv')
colnames(post_tibble) <- c("bu","bl","pu","pl","gamma","D","a")
post_tibble['type'] <- 'Posterior'
post_tibble['type'] <- as.factor(post_tibble$type)

#  Lower and upper boundaries for priors
lm.low<-c(1e-8,1e-7,1e-4,1e-4,1e-6,1e-3,1e-3)
lm.upp<-c(1e-6,1e-5,1,1,2e-3,1,1)

prior_tibble <- tibble(bu = rlunif(1000,min=lm.low[1], max=lm.upp[1],base=10),
                       bl = rlunif(1000,min=lm.low[2], max=lm.upp[2],base=10),
                       pu = rlunif(1000,min=lm.low[3], max=lm.upp[3],base=10),
                       pl = rlunif(1000,min=lm.low[4], max=lm.upp[4],base=10),
                       gamma = rlunif(1000,min=lm.low[5], max=lm.upp[5],base=10),
                       D = rlunif(1000,min=lm.low[6], max=lm.upp[6],base=10),
                       a = rlunif(1000,min=lm.low[7], max=lm.upp[7],base=10),
                )
prior_tibble['type'] <- as.factor('Prior')
post_density <- ggpairs(data=post_tibble,
                        columnLabels = c("beta[U]","beta[L]","p[U]","p[L]","gamma","D","a"), labeller="label_parsed",
                        upper = list(continuous = ggally_density),
                        lower = list(continuous = ggally_density),
                        diag = list(continuous = ggally_densityDiag)
                        ) + theme(panel.spacing = unit(2, "lines"))

test <- bind_rows(post_tibble,prior_tibble)
post_density_prior <- ggpairs(data=test, columns = 1:7, mapping = aes(color = test$type),
                        columnLabels = c("beta[U]","beta[L]","p[U]","p[L]","gamma","D","a"), labeller="label_parsed",
                        upper = list(continuous = ggally_density, rescale = TRUE),
                        lower = list(continuous = ggally_density, rescale = TRUE),
                        diag = list(continuous = ggally_densityDiag, rescale = TRUE)
) + theme(panel.spacing = unit(2, "lines"))

ggsave("Plots/Param_posterior_adaptive_prior.pdf",plot = post_density_prior, width = 40,  height = 40,  units = "cm")
ggsave("Plots/Param_posterior_adapative.pdf",plot = post_density, width = 40,  height = 40,  units = "cm")

early_peak <- read_csv('Results/early_peak_indices.csv',col_names=0)
late_peak <- read_csv('Results/late_peak_indices.csv',col_names=0)
groups <- rep('early',1000)
groups[t(late_peak)] <- 'late'
post_tibble['peak_time'] <- as.factor(groups)
post_density_peak <- ggpairs(data=post_tibble, columns = 1:7, mapping = aes(color = post_tibble$peak_time),
                              columnLabels = c("beta[U]","beta[L]","p[U]","p[L]","gamma","D","a"), labeller="label_parsed",
                              upper = list(continuous = ggally_density, rescale = TRUE),
                              lower = list(continuous = ggally_density, rescale = TRUE),
                              diag = list(continuous = ggally_densityDiag, rescale = TRUE)
) + theme(panel.spacing = unit(2, "lines"))
ggsave("Plots/Param_posterior_adaptive_peak.pdf",plot = post_density_peak, width = 40,  height = 40,  units = "cm")
