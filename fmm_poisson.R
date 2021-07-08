#fmm_poisson.R
## This script implements a Gibbs sampler to estimate a finite mixture of Poisson distributions
## See Fruhwirth-Schnatter, Finite mixture and markov switching models (Chapter 9) for further technical details
## Example uses simulated data

load('data/simulated.data_poisson.rda')

## dependencies: MCMCpack (for random value generation from Dirichlet)
require(MCMCpack)

#Initialisation
K<-2 #number of clusters
n=length(y)
z<-sample(1:K,n,replace=T)
nk<-sapply(1:K,function(k) sum(z==k)) #number of observations assigned to each cluster

#mixture weights (eta)
eta = rdirichlet(1, rep(1,K))

#cluster-specific means (theta)
a0 = 1
b0 = a0/mean(y)

theta = rgamma(K,shape=a0,rate=b0)

#set up MCMC information
MCMC<-list(niter=2000,thin=1,counter=1)
MCMC.traces<-list(z=array(0,c(n,MCMC$niter/MCMC$thin)),eta=array(0,c(K,MCMC$niter/MCMC$thin)),theta=array(0,c(K,MCMC$niter/MCMC$thin)))

#set up progress bar
pbc<-1;pb <- winProgressBar(title="MCMC progress", label="0% done", min=0, max=100, initial=0)
for (t in 1:MCMC$niter){
  
  #Step 1: update latent variable/unobserved memebership to K clusters (z)
  ## note: posterior probabilities of membership calculated on log scale to avoid overflow problems (log-exp-sum trick)
  z <- rep(0,n)
  for (i in 1:n){
    zprob = sapply(1:K,function(k) log(eta[k]) + y[i]*log(theta[k])-theta[k]) #likelihood on log scale. Could also use  dpois(y[i],lambda = theta[k],log=T)
    zprob = exp(zprob-max(zprob)) 
    z[i]<-1+sum(runif(1)>cumsum(zprob/sum(zprob)))
  }
  
  #Step 2: conditional on z, update eta
  nk = sapply(1:K,function(k) sum(z==k))
  eta = rdirichlet(1, rep(1,K)+nk)
  
  #Step 3: conditional on z, update theta
  theta<-sapply(1:K,function(k) rgamma(1,shape=a0+sum(y[z==k]),rate=b0+nk[k]))
  
  #store values
  #store traces
  if ((t/MCMC$thin)==MCMC$counter){
    MCMC.traces$z[,MCMC$counter]<- z
    MCMC.traces$eta[,MCMC$counter] <- eta
    MCMC.traces$theta[,MCMC$counter]<-theta
    MCMC$counter<-MCMC$counter+1
  }
  #update progress bar
  if ((t/pbc)==100){
    Sys.sleep(0.1)
    info <- sprintf("%d%% done", round((t/MCMC$niter)*100))
    setWinProgressBar(pb,100*t/MCMC$niter,label=info)
    pbc<-pbc+1}
}
close(pb)

#save traces
save(MCMC.traces,file='mcmc output/fmm_pois_MCMCoutput.rda')
