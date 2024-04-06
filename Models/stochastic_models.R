# Define a transition function for each transition
# T->E->I->D (upper)
# T->E->I->D (lower)
# Production, clearance, infection, transport, transport
# Production, clearance, infection, immune, transport, transport
# Innate, death
# Adaptive, death
transitions = list(c(Tu = -1, Eu = +1), c(Eu = -1, Iu = +1), c(Iu = -1, Du = +1), c(Du = -1, Tu = +1),
                   c(Tl = -1, El = +1), c(El = -1, Il = +1), c(Il = -1, Dl = +1), c(Dl = -1, Tl = +1),
                   c(Vu = +1), c(Vu = -1), c(Vu = -1), c(Vu = -1), c(Vu = +1),
                   c(Vl = +1), c(Vl = -1), c(Vl = -1), c(Vl = -1), c(Vl = -1), c(Vl = +1),
                   c(Fa = +1), c(Fa = -1),
                   c(X = +1), c(X = -1))
# Define a rate function for:
# T->E->I->D (upper)
# T->E->I->D (lower)
# Production, clearance, infection, transport, transport
# Production, clearance, infection, immune, transport, transport
# Innate, death
# Adaptive, death
lvrates <- function(x, para, t){
  return(c(para$b_u*x['Tu']*x['Vu'], para$g*x['Eu'], para$d*x['Iu'], para$lambda_u*x['Du'], 
           para$b_l*x['Tl']*x['Vl'], para$g*x['El'], para$d*x['Il'], para$lambda_l*x['Dl'],
           para$p*x['Iu'], para$c*x['Vu'], para$g*para$b_u*x['Tu']*x['Vu'], para$t_1*x['Vu'], para$t_2*x['Vl'],
           para$p*x['Il']/(1+para$kappa*x['Fa']), para$c*x['Vl'], para$g*para$b_u*x['Tl']*x['Vl'], para$k*x['Vl']*x['X'], para$t_2*x['Vl'], para$t_1*x['Vu'],
           para$w*x['Vl'], para$delta*x['Fa'],
           para$f*x['Vl'], para$r*x['X']))
}
# Initialise total time and final sizes vector
total_time <- 0
test_cao <- c()
# Set seed for reproducability
set.seed(5536239)
NRuns <- 1
para <- list("lambda_u" = 0, "lambda_l"=0.015, "b_u"=1e-7, "b_l"=1e-7, "g"=4, "c"=2, "d"=10, "p"=1,
             "kappa"=0.045, "gamma"=0.13, "k"=10, "w"=10, "delta"=1, "f"=1e-6, "r"=0.27, "t_1"=1, "t_2"=1)
inits <- c(Tu=1e9, Eu=0, Iu=0, Du=0, Vu=4e4, Tl=1e9, El=0, Il=0, Dl=0, Vl=0, Fa=0, X=0)
maxtime <- 30

library(adaptivetau)
for(i in 1:NRuns){
  # Run current realisation and record final size
  # Also time how long the algorithm takes
  start_time <- Sys.time()
  res = ssa.adaptivetau(inits, transitions, lvrates, para,
                        maxtime, tl.params = list(epsilon=0.05))
  end_time <- Sys.time()
  total_time <- total_time + (end_time-start_time)
  last_ind = dim(res)[1]
}
