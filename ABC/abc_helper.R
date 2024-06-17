library(deSolve)

##############################################################################
#Define the ODE system functions
#############################################################################
diff_host_model <- function(t, pop, para){
  
  #Assign the population matrix into the classes
  Tu <- pop[1]
  Eu <- pop[2]
  Iu <- pop[3]
  Vu <- pop[4]
  Tl <- pop[5]
  El <- pop[6]
  Il <- pop[7]
  Vl <- pop[8]
  X <- pop[9]
  
  #Write down the ODE system
  dTu <- -para$bu*Tu*Vu
  dEu <- para$bu*Tu*Vu - para$g*Eu
  dIu <- para$g*Eu - para$d*Iu
  dVu <- para$pu*Iu - para$c*Vu - para$gamma*para$bu*Tu*Vu - para$D*(Vu-Vl) + para$a*Vl
  dTl <- -para$pl*Tl*Vl
  dEl <- para$bl*Tl*Vl - para$g*El
  dIl <- para$g*El - para$d*Il
  dVl <- para$pl*Il - para$c*Vl - para$gamma*para$bl*Tl*Vl + para$D*(Vu-Vl) - para$a*Vl - para$k*Vl*X
  dX <- para$f*Vl + para$r*X
  
  #Return the derivatives together in a list 
  return(list(c(dTu, dEu, dIu, dVu, dTl, dEl, dIl, dVl, dX)))
  
}

ODE_host_model <- function(para, ICs) {
  
  #Time points for simulation
  t_seq <- seq(0, 10, by = 0.01)
  
  #Solve the ODEs using the ode function from deSolve package
  result <- ode(y = ICs, times = t_seq, func = diff_host_model, parms = para)#, method = "ode45")
  
  #Convert the result to a data frame
  Classes <- data.frame(result)
  
  return(Classes)
}

para <- list("bu" = 4.066e-7, "bl" = 3.67e-7, "g" = 4, "c" = 2, "d" = 5.2, "pu" = 0.298, 
             "pl" = 0.0964, "gamma" = 0.00346, "k" = 20, "f" = 0.56*2.8e-7/7, "r" = 0.27/7, "D"=0.213, "a"=0.147)

##############################################################################
# Define the distance function and a model run function
#############################################################################
dist_fun <- function(theta){
  ICs <- c("Tu" = 4e8,"Eu"=0,"Iu"=0,"Vu"=1.3e3,"Tl"=6.25e9,"El"=0,"Il"=0,"Vl"=0,"X"=0)
  para <- list("bu" = theta[1], "bl" = theta[2], "g" = 4, "c" = 2, "d" = 5.2, "pu" = theta[3], 
               "pl" = theta[4], "gamma" = theta[5], "k" = 20, "f" = 0.56*2.8e-7/7, 
               "r" = 0.27/7, "D"=theta[6], "a"=theta[7])
  ts <- c(401,401,501,501,501,601,601,601,601,601,601,701,701,701,701,801)
  Classes <- ODE_host_model(para, ICs)
  mresults <- Classes$Vu[ts]
  y = 10^c(7.2,5.40,5.1,7.05,7.7,7.68,7.85,7.73,7.01,6.073,4.45,5.894,5.5356,4.325,3.959,6.548)
  return(sum((y-mresults)^2))
}

log_dist_fun <- function(theta){
  ICs <- c("Tu" = 4e8,"Eu"=0,"Iu"=0,"Vu"=1.3e3,"Tl"=6.25e9,"El"=0,"Il"=0,"Vl"=0,"X"=0)
  para <- list("bu" = theta[1], "bl" = theta[2], "g" = 4, "c" = 2, "d" = 5.2, "pu" = theta[3], 
               "pl" = theta[4], "gamma" = theta[5], "k" = 20, "f" = 0.56*2.8e-7/7, 
               "r" = 0.27/7, "D"=theta[6], "a"=theta[7])
  ts <- c(401,401,501,501,501,601,601,601,601,601,601,701,701,701,701,801)
  Classes <- ODE_host_model(para, ICs)
  mresults <- Classes$Vu[ts]
  y = 10^c(7.2,5.40,5.1,7.05,7.7,7.68,7.85,7.73,7.01,6.073,4.45,5.894,5.5356,4.325,3.959,6.548)
  return(sum((log(y)-log(mresults))^2))
}

run_model <- function(theta){
  ICs <- c("Tu" = 4e8,"Eu"=0,"Iu"=0,"Vu"=1.3e3,"Tl"=6.25e9,"El"=0,"Il"=0,"Vl"=0,"X"=0)
  para <- list("bu" = theta[1], "bl" = theta[2], "g" = 4, "c" = 2, "d" = 5.2, "pu" = theta[3], 
               "pl" = theta[4], "gamma" = theta[5], "k" = 20, "f" = 0.56*2.8e-7/7, 
               "r" = 0.27/7, "D"=theta[6], "a"=theta[7])
  ts <- c(401,401,501,501,501,601,601,601,601,601,601,701,701,701,701,801)
  Classes <- ODE_host_model(para, ICs)
  return(Classes$Vu)
}

##############################################################################
#Define other functions
#############################################################################

# Perturbation kernel 
rK <- function(mean, sigma){   
  return(rtmvnorm(1,mean=mean, sigma=sigma, lower=lm.low, upper=lm.upp)) 
}

#  Identity function: H(x)= 1 if x=T
H <- function(x) as.numeric(x>0)

#  Test if prior is non zero
prior.non.zero<-function(par){
  prod(sapply(1:7, function(a) H(par[a]-lm.low[a])* H(lm.upp[a]-par[a])))
}

Norm.Eucl.dist<-function(p1,p2){
  sqrt(sum(((p1-p2)/(lm.upp-lm.low))^2)) }

#  Covariance based on M neighbours
getSigmaNeighbours<-function(M, theta, Theta){
  dist<- sapply(1:N, function(a) Norm.Eucl.dist(as.numeric(theta), as.numeric(Theta[a,])))
  temp<-data.frame(no=seq(1,N), dist)
  temp<-temp[order(temp$dist),]
  sigma<-cov(Theta[temp$no[1:(M+1)],])
  return(sigma)
}

















#Could use Fsolve to find equilibrium
#Define initial conditions need to start close so it doesn't find the non-endemic state
ICs <- c("Tu" = 4e8,"Eu"=0,"Iu"=0,"Vu"=1.3e3,"Tl"=6.25e9,"El"=0,"Il"=0,"Vl"=0,"X"=0)