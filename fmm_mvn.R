#fmm_mvn.R
## This script implements a Gibbs sampler to estimate a finite mixture of Multivariate Normal distributions
## See Fruhwirth-Schnatter, Finite mixture and markov switching models (Chapter 6) for further technical details


## dependencies: MCMCpack, mvtnorm (for random value generation from Dirichlet, MVN, Inverse Wishart distributions)

require(MCMCpack)
require(mvtnorm)


#input data are principal components applied to sampled spikes
## dat: raw data (sampled waveforms)
## dat.pca: processed data following application of PCA
load('data/spikes.rda')

#observed data
y<-as.matrix(dat.pca)

r<-dim(y)[1]
n<-dim(y)[2]

#set hyperparameters; data-driven to help sample behave
hypers<-list(b0=apply(y,1,mean),N0=1,c0=r+1,C0=0.75*cov(t(y)))


#Initialisation
K<-4 #number of clusters
tmp<-kmeans(t(y),K) #centers from k-means as initial 'good guess'
z<-as.vector(tmp$cluster) #latent classification
nk<-sapply(1:K,function(k) sum(z==k)) #number of observations assigned to each cluster

#mixture weights (eta)
eta = rdirichlet(1, rep(1,K))
#cluster means (mu)
mu<-t(tmp$centers)
#cluster var-covar matrices (Sigma)
Sigma<-array(0,c(r,r,K))

for(k in 1:K){Sigma[,,k]<-with(hypers,riwish(c0,C0))}


#set up MCMC information
MCMC<-list(niter=500,thin=1,counter=1)
MCMC.traces<-list(z=array(0,c(n,MCMC$niter/MCMC$thin)),eta=array(0,c(K,MCMC$niter/MCMC$thin)),mu=array(0,c(r,K,MCMC$niter/MCMC$thin)))

#set up progress bar
pbc<-1;pb <- winProgressBar(title="MCMC progress", label="0% done", min=0, max=100, initial=0)
for (t in 1:MCMC$niter){
  
  #Step 1: update latent variable/unobserved memebership to K clusters (z)
  ## note: posterior probabilities of membership calculated on log scale to avoid overflow problems (log-exp-sum trick)
  z <- rep(0,n)
  for (i in 1:n){
    zprob = sapply(1:K, function(k) log(eta[k])+dmvnorm(y[,i], mu[,k], Sigma[,,k], log = TRUE))
    zprob = exp(zprob-max(zprob)) 
    z[i]<-1+sum(runif(1)>cumsum(zprob/sum(zprob)))
  }
  
  #Step 2: conditional on z, update eta
  nk = sapply(1:K,function(k) sum(z==k))
  eta = rdirichlet(1, rep(1,K)+nk)
  
  #Step 3: conditional on z, update mu, sigma
  Sigma<-array(0,c(r,r,K))
  mu<-array(0,c(r,K))
  for (k in 1:K){
    if(nk[k]>1){ybar<-apply(y[,z==k],1,mean)}
    if(nk[k]==1){ybar<-y[,z==k]} #singleton cluster
    mu[,k]<-rmvnorm(1,(nk[k]*ybar+with(hypers,N0*b0))/(nk[k]+hypers$N0),Sigma[,,k]/(hypers$N0+nk[k]))
    D<-as.matrix(y[,z==k]-ybar)
    Sigma[,,k]<-riwish(hypers$c0+nk[k],hypers$C0+with(hypers,(nk[k]*N0/(nk[k]+N0))*(ybar-b0)%*%t(ybar-b0)+D%*%t(D)))
  }
  
  #store values
  #store traces
  if ((t/MCMC$thin)==MCMC$counter){
    MCMC.traces$z[,MCMC$counter]<- z
    MCMC.traces$eta[,MCMC$counter] <- eta
    MCMC.traces$mu[,,MCMC$counter]<-mu
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
save(MCMC.traces,file='fmm_spikesorting_output.rda')
