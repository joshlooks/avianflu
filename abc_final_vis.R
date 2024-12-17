library(GGally)
library(tidyverse)
library(here)
require(latex2exp)
library(KScorrect)
post_tibble <- read_csv('log_adaptive_added_10.csv')
colnames(post_tibble) <- c("bu","bl","pu","pl","gamma","D","a")

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

ggsave("Plots/Param_logposterior_final.pdf",plot = p, width = 40,  height = 40,  units = "cm")

calc_credible_interval <- function(x) {
  quantile(x, probs = c(0.025, 0.975))
}
ci_results <- sapply(post_tibble, calc_credible_interval)
print(ci_df)
