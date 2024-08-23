library(GGally)
library(tidyverse)
library(here)
require(latex2exp)
library(KScorrect)
post_tibble <- read_csv('Posterior/log_adaptive_4.csv')
colnames(post_tibble) <- c("bu","bl","pu","pl","gamma","D","a")
log_post_tibble <- log10(post_tibble)
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
log_prior_tibble <- log10(prior_tibble)
prior_tibble['type'] <- as.factor('Prior')
log_prior_tibble['type'] <- as.factor('Prior')
post_density <- ggpairs(data=post_tibble,
                        columnLabels = c("beta[U]","beta[L]","p[U]","p[L]","gamma","D","a"), labeller="label_parsed",
                        upper = list(continuous = ggally_density),
                        lower = list(continuous = ggally_density),
                        diag = list(continuous = ggally_densityDiag)
                        ) + theme(panel.spacing = unit(2, "lines"))

post_density <- ggpairs(data=post_tibble,
                        columnLabels = c("beta[U]","beta[L]","p[U]","p[L]","gamma","D","a"), labeller="label_parsed",
                        upper = list(continuous = wrap('density',fill = after_stat('level'), geom = "polygon", colour="white")),
                        lower = list(continuous = wrap('density',fill = 'level', geom = "polygon", colour="white")),
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

log_post_density <- ggpairs(data=log_post_tibble,
                            columnLabels = c("beta[U]","beta[L]","p[U]","p[L]","gamma","D","a"), labeller="label_parsed",
                            upper = list(continuous = ggally_density),
                            lower = list(continuous = ggally_density),
                            diag = list(continuous = ggally_densityDiag)
) + theme(panel.spacing = unit(2, "lines"))
ggsave("Plots/Param_logposterior_adaptive.pdf",plot = log_post_density, width = 40,  height = 40,  units = "cm")

# Function to create contour plots
kde_with_contours <- function(data, mapping, ...){
  ggplot(data = data, mapping = mapping) +
    geom_density2d(aes(color = ..level..)) +
    scale_color_viridis_c()
}

# Create the pair plot with KDE and log-scaled axes
p <- ggpairs(post_tibble, columnLabels = c("beta[U]","beta[L]","p[U]","p[L]","gamma","D","a"), labeller="label_parsed",
             upper = list(continuous = wrap(kde_with_contours)),
             lower = list(continuous = wrap(kde_with_contours)),
             diag = list(continuous = wrap("densityDiag"))) +
  theme(panel.spacing = unit(2, "lines"))

scalelog10<-function(x=7,g){
  # Bottom triangle
  for (i in 2:x){ 
    for (j in 1:(i-1)) { 
      g[i,j]<-g[i,j] + scale_x_continuous(trans='log10') +
        scale_y_continuous(trans='log10')
    } } 
  # Upper triangle
  for (i in 1:(x-1)){ 
    for (j in (i+1):x) { 
      g[i,j]<-g[i,j] + scale_x_continuous(trans='log10') +
        scale_y_continuous(trans='log10')
    } } 
  for (i in 1:x){ #for the diagonal
    g[i,i]<-g[i,i]+ scale_x_continuous(trans='log10')  } 
  return(g) }

p <- scalelog10(x=7,p)

log_post_tibble['type'] <- 'Posterior'
log_post_tibble['type'] <- as.factor(log_post_tibble$type)
log_test <- bind_rows(log_post_tibble,log_prior_tibble)
log_post_density_prior <- ggpairs(data=log_test, columns = 1:7, mapping = aes(color = log_test$type),
                              columnLabels = c("beta[U]","beta[L]","p[U]","p[L]","gamma","D","a"), labeller="label_parsed",
                              upper = list(continuous = ggally_density, rescale = TRUE),
                              lower = list(continuous = ggally_density, rescale = TRUE),
                              diag = list(continuous = ggally_densityDiag, rescale = TRUE)
) + theme(panel.spacing = unit(2, "lines"))

ggsave("Plots/Param_logposterior_adaptive.pdf",plot = p, width = 40,  height = 40,  units = "cm")

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
