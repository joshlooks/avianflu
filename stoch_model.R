# Load libraries
require(here)
require(adaptivetau)
require(matrixStats)
require(tidyverse)
require(latex2exp)

# Load stochastic model
source(paste(here(), "/Models/stochastic_models.R", sep=""))

NRuns <- 1000
TU <- 4e8
TL <- 6.25e9
para <- list("b_u"=1.5e-8, "b_l"=1.5e-6, "g"=4, "c"=2, "d"=5.2, "pu"=5e7/(TU+TL), "pl"=5e7/(TU+TL),
             "gamma"=0, "k"=20, "f"=0.56*2.8e-6, "r"=0.27, "D"=0.01, "a"=0.1)
inits <- c(Tu=TU, Eu=0, Iu=0, Vu=1.3e3, Tl=6.25e9, El=0, Il=0, Vl=0, X=0, Du=0, Dl=0)
maxtime <- 10
times <- seq(0,10,0.05)
Vlres <- matrix(0,nrow=NRuns,ncol=length(times))
Dlres <- matrix(0,nrow=NRuns,ncol=length(times))
set.seed(1234)
for(i in 1:NRuns){
  # Run current realisation
  res = ssa.adaptivetau(inits, transitions_diffadv, lvrates_diffadv, para,
                        maxtime, tl.params = list(epsilon=0.05))
  Vlres[i,] <- approx(res[,1],res[,9], xout=times)$y
  Dlres[i,] <- approx(res[,1],res[,12], xout=times)$y
}

Vlresquants <- colQuantiles(Vlres, probs = c(0.025,0.975))
Vlresmed <- colQuantiles(Vlres, probs = 0.5)
Dlresquants <- colQuantiles(Dlres, probs = c(0.025,0.975))
Dlresmed <- colQuantiles(Dlres, probs = 0.5)

Vltibble <- tibble(times,Vlresquants[,1],Vlresquants[,2],Vlresmed)
colnames(Vltibble) <- c("t","LI","UI","Med")
vlplot <- ggplot(Vltibble, aes(t, Med)) +                                      
  geom_line(size = 1) + 
  geom_ribbon(aes(ymin=LI, ymax=UI), alpha=0.1, fill = rgb(0.502,0,0,0.5),  
              color = rgb(0.502,0,0), linetype = "dotted") + 
  labs(x = "Time (days)",
       y = TeX("Viral load (TCID$_{50}$)/ml")) +
  scale_y_continuous(trans='log10') + theme_bw(base_size=20)
ggsave("Plots/Vlstoch.pdf",plot = vlplot)

Dltibble <- tibble(times,Dlresquants[,1],Dlresquants[,2],Dlresmed)
colnames(Dltibble) <- c("t","LI","UI","Med")
dlplot <- ggplot(Dltibble, aes(x=t,y=Med)) +                                      
  geom_line(size = 1) + 
  geom_ribbon(aes(ymin=LI, ymax=UI), alpha=0.1, fill = rgb(0.502,0,0,0.5),
              color = rgb(0.502,0,0), linetype = "dotted") +
  scale_y_continuous(trans='log10') + theme_bw(base_size=20) + 
  labs(x = "Time (days)",
       y = TeX("Dead cells"))
ggsave("Plots/Dlstoch.pdf",plot = dlplot)
