library(GGally)
library(tidyverse)
library(here)
source(paste0(here(),'/ABC/abc_helper.R'))
post_tibble <- read_csv('test_7.csv')
colnames(post_tibble) <- c("bu","bl","pu","pl","gamma","D","a")

post_density <- ggpairs(data=post_tibble, 
        upper = list(continuous = "density"),
        lower = list(continuous = "density"))
NRuns <- dim(post_tibble)[1]
ts <- c(401,401,501,501,501,601,601,601,601,601,601,701,701,701,701,801)
Vlres <- matrix(0,nrow=NRuns,ncol=length(ts))
for(i in 1:NRuns){
  # Run current realisation
  theta = post_tibble[i,]
  res = run_model(theta)
  ts <- c(401,401,501,501,501,601,601,601,601,601,601,701,701,701,701,801)
  Vlres[i,] <- res$Vu[ts]
}

data_days <- c(4,4,5,5,5,6,6,6,6,6,6,7,7,7,7,8)
days_loads <- 10^c(7.2,5.40,5.1,7.05,7.7,7.68,7.85,7.73,7.01,6.073,4.45,5.894,5.5356,4.325,3.959,6.548)
Datatibble <- tibble(data_days, days_loads)
colnames(Datatibble) <- c("t","load")

Vuresquants <- colQuantiles(Vures, probs = c(0.025,0.975))
Vuresmed <- colQuantiles(Vures, probs = 0.5)

Vutibble <- tibble(times,Vuresquants[,1],Vuresquants[,2],Vuresmed)
colnames(Vutibble) <- c("t","LI","UI","Med")
vuplot <- ggplot(Vutibble, aes(t, Med)) +                                      
  geom_line(size = 1) + 
  geom_ribbon(aes(ymin=LI, ymax=UI), alpha=0.1, fill = rgb(0.502,0,0,0.5),  
              color = rgb(0.502,0,0), linetype = "dotted") + 
  labs(x = "Time (days)",
       y = TeX("Viral load (TCID$_{50}$)/ml")) +
  scale_y_continuous(trans='log10') + theme_bw(base_size=20) +
  geom_point(aes(x=data_days,y=days_loads),color="#6680FF",data=Datatibble)
ggsave("Plots/Vu_post.pdf",plot = vlplot, width = 20,  height = 10,  units = "cm")