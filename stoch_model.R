# Load libraries
library(here)
library(adaptivetau)

# Load stochastic model
source(paste(here(), "/Models/stochastic_models.R", sep=""))

# Initialise total time and final sizes vector
total_time <- 0
test_cao <- c()
NRuns <- 20
para <- list("lambda_u" = 0.015, "lambda_l"=0.015, "b_u"=1.9e-9, "b_l"=1.9e-7, "g"=4, "c"=2, "d"=10, "pu"=1, "pl"=1,
             "kappa"=0.045, "gamma"=0.13, "k"=10, "w"=10, "delta"=1, "f"=1e-6, "r"=0.27, "D"=1, "a"=0.01)
inits <- c(Tu=1e9, Eu=0, Iu=0, Vu=4e4, Tl=1e9, El=0, Il=0, Vl=0, X=0, Du=0, Dl=0)
maxtime <- 30
maxVl <- rep(0,NRuns)
maxVu <- rep(0,NRuns)
Du <- rep(0,NRuns)
Dl <- rep(0,NRuns)
for(i in 1:NRuns){
  # Run current realisation and record final size
  # Also time how long the algorithm takes
  start_time <- Sys.time()
  res = ssa.adaptivetau(inits, transitions_diffadv, lvrates_diffadv, para,
                        maxtime, tl.params = list(epsilon=0.05))
  maxVl[i] <- max(res[,9])
  maxVu[i] <- max(res[,5])
  Du[i] <- res[dim(res)[1],11]
  Dl[i] <-res[dim(res)[1],12]
}