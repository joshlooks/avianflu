library(tmvtnorm)
library(here)
library(KScorrect)
source(paste(here(),'/ABC/abc_helper.R',sep=""))

# Number of particles
N <- 1000

# Number of neighbours for covariance matrix calculations
M <- 50

# Epsilon values for temporal data 
# etol <- c(1e18,5e17,1e17,5e16,2e16,1.5e16,1.35e16)
etol <- c(150,145,140,135,130,127.5,125)  

# Number of generations
G <- length(etol)

# Number of simulations for each parameter set
n <- 1

#  Lower and upper boundaries for priors
lm.low<-c(1e-8,1e-7,1e-4,1e-4,1e-6,1e-3,1e-3)
lm.upp<-c(1e-6,1e-5,1,1,2e-3,1,1)


# Empty matrices to store results (population plus 5 model parameters)
res.old<-matrix(ncol=7,nrow=N)
res.new<-matrix(ncol=7,nrow=N)

# Empty vectors to store weights
w.old<-matrix(ncol=1,nrow=N)
w.new<-matrix(ncol=1,nrow=N)

for(g in 1:G){  
  
  #Initiate counter
  i<-1	
  print(paste0('Generation: ', g, "beginning"))
  while(i <= N){ # While the number of accepted particles is less than N_particles
    if(g==1){
      # Sample from prior distributions 
      bu_star<- rlunif(1,min=lm.low[1], max=lm.upp[1],base=10)
      bl_star<- rlunif(1,min=lm.low[2], max=lm.upp[2],base=10)
      pu_star<- rlunif(1,min=lm.low[3], max=lm.upp[3],base=10)
      pl_star<- rlunif(1,min=lm.low[4], max=lm.upp[4],base=10)
      gamma_star<- rlunif(1,min=lm.low[5], max=lm.upp[5],base=10)
      D_star<- rlunif(1,min=lm.low[6], max=lm.upp[6],base=10)
      a_star<- rlunif(1,min=lm.low[7], max=lm.upp[7],base=10)
      par <- c(bu_star, bl_star, pu_star, pl_star, gamma_star, D_star, a_star)
    } else {
      #  Select particle from previous generation
      p<-sample(seq(1,N),1,prob=w.old)		
      sigma<-Sigma[[p]]
      par<- rK(res.old[p,],sigma)
    }
    #  Test if prior non zero
    if(prior.non.zero(par)) {
      # Set number of accepted simulations to zero
      m<-0
      # D_star <- dist_fun(par)
      D_star <- log_dist_fun(par)
      if (D_star <= etol[g]){
        # Store results
        res.new[i,]<-par  
        # Calculate weights
        w1<-prod(sapply(1:7, function(b) dlunif(res.new[i,b], min=lm.low[b], max=lm.upp[b],base=10)))
        if(g==1){
          w2<-1
        } else {
          w2<-sum(sapply(1:N, function(a) w.old[a]* dtmvnorm(res.new[i,], mean=res.old[a,], sigma=sigma, lower=lm.low, upper=lm.upp)))
        }
        w.new[i] <- w1/w2
        # Update counter
        i <- i+1
        if (i %% 100 == 0){
          print(paste0('Particle ', i, " done"))
        }
      }
    } 
  }
  Sigma <- list(NA, N)
  for(p in 1:N){
    Sigma[[p]]<- getSigmaNeighbours(M, res.new[p,], res.new) 
  }
  res.old<-res.new
  w.old<-w.new/sum(w.new)
  
  write.csv(res.new, file = paste("log_",g,".csv",sep=""), row.names=FALSE)
}
