#fmm_BUGS.R
# Example of univariate Normal mixture model using R2OpenBUGS
library(R2OpenBUGS)
data("fish", package = "bayesmix") #fishery data from bayesmix

#define number of mixture components/clusters
K = 2
N = nrow(fish)

b0 = mean(fish$fish)
c0 = 2.5
C0=0.5*var(fish$fish)

dat = list(y = fish$fish, N = N, K=K,alpha = rep(1,K),b0=b0,N0=1,c0=c0,C0=C0)

#plot the data
hist(dat$y,breaks=20)

#initialise classifications. For BUGS, treat as missing but assign at least one observation to each cluster as initial values

z = rep(NA,N)
z[1] = 1
#z[N] = 2
dat$z = z

#set initial values
theta.init = kmeans(dat$y,K)$centers[,1]
inits <- function(){list(theta = theta.init, tau = rep(1,K))}

#set parameters to track
parameters <- c("theta", "eta", "sigma","z")

fmm.sim <- bugs(data=dat, inits=inits, parameters=parameters, model.file = 'fmmBUGS.txt',
                n.chains = 1,n.burnin = 10000, n.iter = 20000,
                working.directory = NULL,debug=F)

#review the output
eta.mcmc = data.frame(fmm.sim$sims.list$eta)
theta.mcmc = data.frame(fmm.sim$sims.list$theta)

