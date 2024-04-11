# Define a transition function for each transition
# T->E->I-> (upper)
# T->E->I-> (lower)
# Production, clearance, infection, transport, transport
# Production, clearance, infection, immune
# Adaptive, growth
transitions = list(c(Tu = -1, Eu = +1), #Tu -> Eu
                   c(Eu = -1, Iu = +1), #Eu -> Iu
                   c(Iu = -1),          #Iu -> D
                   c(Tl = -1, El = +1), #Tl -> El
                   c(El = -1, Il = +1), #El -> Il
                   c(Il = -1),          #Il -> D
                   c(Vu = +1),          #Production from Iu
                   c(Vu = -1),          #Clearance of Vu
                   c(Vu = -1),          #Infection from Vu
                   c(Vu = -1, Vl = +1), #Vu -> Vl transport
                   c(Vu = +1, Vl = -1), #Vl -> Vu transport
                   c(Vl = +1),          #Production from Il
                   c(Vl = -1),          #Clearance of Vl
                   c(Vl = -1),          #Infection from Vl
                   c(Vl = -1),          #Immune response to Vl
                   c(X = +1),           #Adaptive growth from Vl
                   c(X = +1))           #Adaptive exponential growth
# Define a rate function for:
# T->E->I (upper)
# T->E->I (lower)
# Production, clearance, infection, transport, transport
# Production, clearance, infection, immune
# Adaptive, growth
lvrates <- function(x, para, t){
  return(c(infu=para$b_u*x['Tu']*x['Vu'],             #Tu -> Eu
           onu=para$g*x['Eu'],                        #Eu -> Iu
           Du=para$d*x['Iu'],                         #Iu -> D
           infl=para$b_l*x['Tl']*x['Vl'],             #Tl -> El
           onl=para$g*x['El'],                        #El -> Il
           Dl=para$d*x['Il'],                         #Il -> D
           pru=para$pu*x['Iu'],                       #Production from Iu
           clu=para$c*x['Vu'],                        #Clearance of Vu
           infvu=para$gamma*para$b_u*x['Tu']*x['Vu'], #Infection from Vu
           tu=para$t_1*x['Vu'],                       #Vu -> Vl transport
           tl=para$t_2*x['Vl'],                       #Vl -> Vu transport
           prl=para$pl*x['Il'],                       #Production from Il
           cll=para$c*x['Vl'],                        #Clearance of Vl
           infvl=para$gamma*para$b_l*x['Tl']*x['Vl'], #Infection from Vl
           iml=para$k*x['Vl']*x['X'],                 #Immune response to Vl
           adx=para$f*x['Vl'],                        #Adaptive growth from Vl
           gwx=para$r*x['X']))                        #Adaptive exponential growth
}
# Initialise total time and final sizes vector
total_time <- 0
test_cao <- c()
NRuns <- 30
para <- list("lambda_u" = 0, "lambda_l"=0.015, "b_u"=1e-7, "b_l"=1e-7, "g"=4, "c"=10, "d"=2, "pu"=1, "pl"=1,
             "kappa"=0.045, "gamma"=0.13, "k"=10, "w"=10, "delta"=1, "f"=1e-6, "r"=0.27, "t_1"=1, "t_2"=1)
inits <- c(Tu=1e9, Eu=0, Iu=0, Vu=4e4, Tl=1e9, El=0, Il=0, Vl=0, X=0)
maxtime <- 30
maxVl <- rep(0,NRuns)
maxVu <- rep(0,NRuns)
library(adaptivetau)
for(i in 1:NRuns){
  # Run current realisation and record final size
  # Also time how long the algorithm takes
  start_time <- Sys.time()
  res = ssa.adaptivetau(inits, transitions, lvrates, para,
                        maxtime, tl.params = list(epsilon=0.05))
  maxVl[i] <- max(res[,9])
  maxVu[i] <- max(res[,5])
}
