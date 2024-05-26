library(GGally)
library(tidyverse)
library(here)
require(latex2exp)
source(paste0(here(),'/ABC/abc_helper.R'))
post_tibble <- read_csv('Posterior/log_7.csv')
colnames(post_tibble) <- c("bu","bl","pu","pl","gamma","D","a")
post_density <- ggpairs(data=post_tibble,
                        columnLabels = c("beta[U]","beta[L]","p[U]","p[L]","gamma","D","a"), labeller="label_parsed",
                        upper = list(continuous = "density"),
                        lower = list(continuous = "density"),
                        ) + theme(panel.spacing = unit(2, "lines"))
ggsave("Plots/Param_posterior.pdf",plot = post_density, width = 40,  height = 40,  units = "cm")
