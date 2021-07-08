#slice sampler for DP mixture model.
#Nicole White, 2015
#email: nm.white@qut.edu.au


#This script allows for the estimation of a Multivariate Gaussian DPM using Markov chain Monte Carlo (MCMC).

#Implementation requires the following input/s:
#  y : matrix of observations (PCs corresponding to each sampled spike).  
# Dimension is r (number of PCs/dimensions) x n (number of spikes/observations).  

#Values for each hyperparameter are specific as in the main text of the paper, in the list 'hypers'
#For an initial number of clusters (K), the latent classification variable (z) is initialised using K-means
#Information on the number of MCMC samples taken is contained in the list 'MCMC'.  Corresponding samples drawn at each (thinned) iteration are stored in 'MCMC.traces'.

## dependencies: MCMCpack, mvtnorm (for random value generation from Dirichlet, MVN, Inverse Wishart distributions); mcclust for posterior inference

require(MCMCpack)
require(mvtnorm)

#load DPM functions - additional functions to summarise each MCMC step
source("slice sampler DPM functions.R")

load('data/spikes.rda')

#observed data
y<-as.matrix(dat.pca)

r<-dim(y)[1]
n<-dim(y)[2]


#set hyperparameters
hypers<-list(b0=apply(y,1,mean),N0=0.01,c0=r+1,C0=0.75*cov(t(y)),alpha.a=1,alpha.b=1)

#Initialisation
#clustering structure 
K<-10;Kstar<-K
tmp<-kmeans(t(y),K)
z<-as.vector(tmp$cluster)
nk<-sapply(1:K,function(k) sum(z==k))

#cluster means (mu)
mu<-t(tmp$centers)

#cluster var-covar matrices (Sigma)
Sigma<-array(0,c(r,r,K))
for(k in 1:K){Sigma[,,k]<-with(hypers,riwish(c0,C0))}

#DP concentration parameter (alpha)
alpha<-1

#stick breaking weights (v) and cluster weights (weights)
v<-rep(0,K);weights<-v
v[1]<-rbeta(1,1,alpha)
weights[1]<-v[1]
for(k in 2:K){
v[k]<-rbeta(1,1,alpha)
weights[k]<-v[k]*prod(1-v[1:(k-1)])
}
rm(tmp)



#set up MCMC information
MCMC<-list(niter=1000,thin=1,counter=1)
MCMC.traces<-list(z=array(0,c(n,MCMC$niter/MCMC$thin)),alpha=rep(0,MCMC$niter/MCMC$thin),K=rep(0,MCMC$niter/MCMC$thin))

#set up progress bar
pbc<-1;pb <- winProgressBar(title="MCMC progress", label="0% done", min=0, max=100, initial=0)
for (t in 1:MCMC$niter){
#STEP 1: Determine current number of alive components
Kstar<-max(z);nk<-nk[1:Kstar]

#STEP 2: Update stick breaking weights for alive components
v<-vupdate(Kstar,nk,alpha,n)
weights<-weightsupdate(Kstar,v)

#STEP 3: Update MVN parameters for alive components
out<-muSigmaupdate(y,z,nk,hypers,r,Kstar)
mu<-out$mu
Sigma<-out$Sigma

#STEP 4:label switching move
out<-labelswitch(Kstar,z,nk,weights,v,alpha,mu,Sigma)
z<-out$z
nk<-out$nk
v<-out$v
weights<-out$weights
mu<-out$mu
Sigma<-out$Sigma

#STEP 5:Update u
u<-sapply(1:n,function(t) runif(1,0,weights[z[t]]))

#STEP 6: Update alpha
alpha<-updatealpha(n,alpha,Kstar,hypers$alpha.a,hypers$alpha.b)

#STEP 7: Determine if more components are required (potentials)
out<-determineKP(min(u),weights,Kstar,v,alpha)
KP<-out$KP
Ktotal<-Kstar+KP
if(KP>0){v<-out$v; weights<-out$weights}

#STEP 8: Update theta for any potential components
if (KP>0){
mu<-cbind(mu,array(0,c(r,KP)))
for (k in (Kstar+1):Ktotal){
Sigma<-array(c(Sigma,matrix(0,r,r)),c(r,r,k))
Sigma[,,k]<-riwish(hypers$c0,hypers$C0)
mu[,k]<-rmvnorm(1,hypers$b0,Sigma[,,k]/hypers$N0)
} #end of for loop
}#end of if statement

#STEP 9: Update z
out<-zupdate(Ktotal,n,y,u,weights,mu,Sigma)
nk<-out$nk
z<-out$z

#CLEAN UP: order by decreasing z
ztmp<-z;counter<-1;
tmp<-unique(ztmp)
for(i in tmp){
z[ztmp==i]<-counter
counter<-counter+1
}
weights<-weights[tmp]
v<-v[tmp]
mu<-mu[,tmp]
Sigma<-Sigma[,,tmp]
nk<-nk[tmp]
rm(tmp)

#store traces
if ((t/MCMC$thin)==MCMC$counter){
MCMC.traces$z[,t]<-norm.label(z) #relabelled to 1:length(unique(z))
MCMC.traces$alpha[t]<-alpha
MCMC.traces$K[t]<-sum(nk>0) #number of occupied components
MCMC$counter<-MCMC$counter+1
}

#update progress bar
if ((t/pbc)==10){
Sys.sleep(0.1)
info <- sprintf("%d%% done", round((t/MCMC$niter)*100))
setWinProgressBar(pb,t/MCMC$niter,label=info)
pbc<-pbc+1}

} #end of t loop
close(pb)

#Post-processing: Find MAP clustering based on PEAR.
sim.mean<-comp.psm(t(MCMC.traces$z[,1:t-1]))
image(1:n,1:n,sim.mean[order(y[1,]),order(y[1,])])
zmap<-maxpear(sim.mean,t(MCMC.traces$z[,1:t-1]),method="draws")$cl 