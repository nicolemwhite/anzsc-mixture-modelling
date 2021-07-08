#fmm_multinomial.R
## This script implements a Gibbs sampler to estimate a finite mixture of Multinomial distributions
## See Fruhwirth-Schnatter, Finite mixture and markov switching models (Chapter 6) for further technical details


## dependencies: MCMCpack (for random value generation from Dirichlet)
require(MCMCpack)

#for this example, we'll work with a simulated dataset
load('data/simulated.data_multinomial.rda')

#observed data (y)

n = nrow(y) #number of subjects
J = ncol(y) #number of items
R = 3 #three levels per item

#Initialisation
K<-3 #number of clusters
z<-sample(1:K,n,replace=T)
nk<-sapply(1:K,function(k) sum(z==k)) #number of observations assigned to each cluster

#mixture weights (eta)
eta = rdirichlet(1, rep(1,K))

#cluster probabilities (theta)
theta = array(0,c(R,J,K))
for (k in 1:K){
  for (j in 1:J){
    theta[,j,k] = rdirichlet(1, rep(1/R,R))
  }
}


#set up MCMC information
MCMC<-list(niter=2000,thin=1,counter=1)
MCMC.traces<-list(z=array(0,c(n,MCMC$niter/MCMC$thin)),eta=array(0,c(K,MCMC$niter/MCMC$thin)),theta=array(0,c(R,J,K,MCMC$niter/MCMC$thin)))

#set up progress bar
pbc<-1;pb <- winProgressBar(title="MCMC progress", label="0% done", min=0, max=100, initial=0)
for (t in 1:MCMC$niter){
  
  #Step 1: update latent variable/unobserved memebership to K clusters (z)
  ## note: posterior probabilities of membership calculated on log scale to avoid overflow problems (log-exp-sum trick)
  z <- rep(0,n)
  for (i in 1:n){
    ymat = sapply(1:R,function(x) as.numeric(y[i,]==x)) #changes multinomial data to binary matrix to save extra loops. Dimension is number items X number of levels (assumed equal here)
    zprob = sapply(1:K,function(k) log(eta[k]) + sum(diag(ymat%*%log(theta[,,k])))) #likelihood on log scale. Could also use dmultinom here; e.g. sum(sapply(1:J,function(j) dmultinom(ymat[j,],prob=theta[,j,k],log=T)))
    zprob = exp(zprob-max(zprob)) 
    z[i]<-1+sum(runif(1)>cumsum(zprob/sum(zprob)))

  }
  
  #Step 2: conditional on z, update eta
  nk = sapply(1:K,function(k) sum(z==k))
  eta = rdirichlet(1, rep(1,K)+nk)
  
  #Step 3: conditional on z, update theta
  theta<-array(0,c(R,J,K))
  for (k in 1:K){
      for (j in 1:J){
        m = sapply(1:R,function(x) sum(y[z==k,j]==x)) #number of observations per level of multinomial (item j, cluster k)
        theta[,j,k] = rdirichlet(1,m+rep(1/R,R))
      }
  }
  
  #store values
  #store traces
  if ((t/MCMC$thin)==MCMC$counter){
    MCMC.traces$z[,MCMC$counter]<- z
    MCMC.traces$eta[,MCMC$counter] <- eta
    MCMC.traces$theta[,,,MCMC$counter]<-theta
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
save(MCMC.traces,file='mcmc output/fmm_multinomal_MCMCoutput.rda')
