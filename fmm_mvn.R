# Multivariate normal mixture, gibbs sampler applied to spike sorting dataset
## See Fruhwieth-Schnatter, Ch6 for further technical details.
require(mcclust)
require(MCMCpack)
require(mvtnorm)
require(tidyverse)

#input data are principal components applied to sampled spikes
load('data/spikes.rda')

y<-as.matrix(dat.pca)

r<-dim(y)[1]
n<-dim(y)[2]

#set hyperparameters
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
MCMC<-list(niter=100,thin=1,counter=1)
MCMC.traces<-list(z=array(0,c(n,MCMC$niter/MCMC$thin)),eta=array(0,c(K,MCMC$niter/MCMC$thin)),mu=array(0,c(r,K,MCMC$niter/MCMC$thin)))

#set up progress bar
pbc<-1;pb <- winProgressBar(title="MCMC progress", label="0% done", min=0, max=100, initial=0)
for (t in 1:MCMC$niter){
  
  #Step 1: update latent variable z
  ## note: posterior probabilities of membership calculated on log scale to avoie overflow problems (log-exp-sum trick)
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

#save traces
save(MCMC.traces,file='fmm_spikesorting_output.rda')


#label switching
##apply prior constraint on eta: only accept draws where eta1<...<etaK
## e.g mu[1,1:K]
eta.post <- mu.post <- rep(0,K)
for (t in 1:MCMC$niter){
  eta_order = order(MCMC.traces$eta[,t])
  if (sum(eta_order==1:K)==K){
    eta.post = rbind(eta.post,MCMC.traces$eta[,t])
    mu.post = rbind(mu.post,MCMC.traces$mu[1,,t])
  }
}

## apply Stephens' relabelling algorithm
z_unswitched = relabel(t(MCMC.traces$z))
z.post = z_unswitched$cls 
zprob.post = z_unswitched$P

## apply similarity matrix to visualise clustering
z.sim = comp.psm(t(MCMC.traces$z))

#reorder for heatmap
tmp<-hclust(dist(z.sim))
z.sim_reorder<-z.sim[tmp$order,tmp$order]
z.sim_long <- reshape2::melt(z.sim_reorder,varnames = c('z1','z2'))
ggplot(z.sim_long,aes(x=z1,y=z2,fill=value))+
  geom_raster()+xlab('')+ylab('')+
  scale_fill_continuous(low='#CCCCCC',high='#333333')+
  scale_x_discrete(expand = c(0, 0)) +scale_y_discrete(expand = c(0, 0)) +
  theme_grey()+theme(legend.title = element_blank(),legend.position='none',plot.title = element_text(hjust = 0.5,size=14),text=element_text(size=14))


image(1:n,1:n,z.sim[order(y[1,]),order(y[1,])],col=grey.colors(25))

