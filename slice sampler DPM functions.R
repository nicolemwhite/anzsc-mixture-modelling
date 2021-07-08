#functions to implement dpm using slice sampler
require(mcclust)
require(MCMCpack)
require(mvtnorm)

determineKP<-function(min.u,weights,Kstar,v,alpha){
counter<-Kstar
if (sum(weights)<(1-min.u)){
while(sum(weights)<(1-min.u)){
v<-c(v,0)
weights<-c(weights,0)
v[counter+1]<-rbeta(1,1,alpha)
weights[counter+1]<-v[counter+1]*prod(1-v[1:counter])
counter<-counter+1
}
}
#number of potential components
KP<-counter-Kstar
output<-list(v=v,weights=weights,KP=KP)
return(output)
}

zupdate<-function(Kstar,n,y,u,weights,mu,Sigma){
nk<-rep(0,Kstar)
for (i in 1:n){
zprob<-rep(0,Kstar)
#zprob[weights>u[i]]<-sapply(which(weights>u[i]),function(k) -0.5*sum(log(eigen(Sigma[,,k])$values))-0.5*(y[,i]-mu[,k])%*%(solve(Sigma[,,k],y[,i]-mu[,k])))
zprob[weights>u[i]]<-sapply(which(weights>u[i]),function(k) dmvnorm(y[,i], mu[,k], Sigma[,,k], log = TRUE))

zprob[weights>u[i]]<-exp(zprob[weights>u[i]]-max(zprob[weights>u[i]]))
z[i]<-1+sum(runif(1)>cumsum(zprob/sum(zprob)))
nk[z[i]]<-nk[z[i]]+1
}
output<-list(nk=nk,z=z)
return(output)
}

vupdate<-function(Kstar,nk,alpha,n){
v<-rep(0,Kstar);weights<-rep(0,Kstar)
v[1]<-rbeta(1,1+sum(z==1),alpha+sum(z>1))
for (k in 2:Kstar){
	v[k]<-rbeta(1,1+nk[k],alpha+n-sum(nk[1:k]))
}
return(v)
}
weightsupdate<-function(Kstar,v){
weights<-rep(0,Kstar)
weights[1]<-v[1]
for (k in 2:Kstar){
	weights[k]<-v[k]*prod(1-v[1:k-1])
}
return(weights)
}

muSigmaupdate<-function(y,z,nk,hypers,P,Kstar){
Sigma<-array(0,c(P,P,Kstar))
mu<-array(0,c(P,Kstar))
for (k in 1:Kstar){
if(nk[k]>1){ybar<-apply(y[,z==k],1,mean)}
if(nk[k]==1){ybar<-y[,z==k]}
D<-as.matrix(y[,z==k]-ybar)
Sigma[,,k]<-riwish(hypers$c0+nk[k],hypers$C0+with(hypers,(nk[k]*N0/(nk[k]+N0))*(ybar-b0)%*%t(ybar-b0)+D%*%t(D)))
mu[,k]<-rmvnorm(1,(nk[k]*ybar+with(hypers,N0*b0))/(nk[k]+hypers$N0),Sigma[,,k]/(hypers$N0+nk[k]))
}

output<-list(mu=mu,Sigma=Sigma)
return(output)
}


#label switching
labelswitch<-function(Kstar,z,nk,weights,v,alpha,mu,Sigma){
move<-runif(1)
if (move<=(1/2)){
#choose two 'alive' components
kchoose<-sample(1:Kstar, 2)
logr<-(nk[kchoose[1]]-nk[kchoose[2]])*(log(weights[kchoose[2]])-log(weights[kchoose[1]]))
}
if(move>(1/2)){
kchoose<-rep(0,2)
#choose any two components from 1:Kstar-1 (ie. empties ok for this move)
kchoose[1]<-sample(1:(Kstar-1),1)
kchoose[2]<-kchoose[1]+1
logr<-nk[kchoose[1]]*log(1-v[kchoose[2]])-nk[kchoose[2]]*log(1-v[kchoose[1]])
}
#if(move>(2/3)){
#tmp<-sample(2:sum(nk>0), 1)
#kchoose<-rep(0,2)
#kchoose[1]<-which(nk>0)[tmp]
#kchoose[2]<-which(nk>0)[tmp-1]
#weightscomb<-weights[kchoose[1]]+weights[kchoose[2]]
#R1<-(1+alpha+nk[kchoose[2]]+sum(nk[(1:Kstar)>kchoose[2]]))/(alpha+nk[kchoose[2]]+sum(nk[(1:Kstar)>kchoose[2]]))
#R2<-(alpha+nk[kchoose[1]]+sum(nk[(1:Kstar)>kchoose[2]]))/(1+alpha+nk[kchoose[1]]+sum(nk[(1:Kstar)>kchoose[2]]))
#weightsdash<-weights[kchoose[2]]*R1+weights[kchoose[1]]*R2
#logr<-(nk[kchoose[1]]+nk[kchoose[2]])*(log(weightscomb/weightsdash))+nk[kchoose[2]]*log(R1)+nk[kchoose[1]]*log(R2)
#}

if (log(runif(1))<min(0,logr)){
ztmp<-z;vtmp<-v;nktmp<-nk;mutmp<-mu;sigmatmp<-Sigma
z[ztmp==kchoose[1]]<-kchoose[2]
z[ztmp==kchoose[2]]<-kchoose[1]
nk[kchoose[1]]<-nktmp[kchoose[2]]
nk[kchoose[2]]<-nktmp[kchoose[1]]
mu[,kchoose[1]]<-mutmp[,kchoose[2]]
mu[,kchoose[2]]<-mutmp[,kchoose[1]]
Sigma[,,kchoose[1]]<-sigmatmp[,,kchoose[2]]
Sigma[,,kchoose[2]]<-sigmatmp[,,kchoose[1]]
if (move>(1/2)){
v[kchoose[1]]<-vtmp[kchoose[2]]
v[kchoose[2]]<-vtmp[kchoose[1]]
weights<-weightsupdate(Kstar,v)
}
#if (move>(2/3)){
#weightstmp<-weights;vtmp<-v
#weights[kchoose[1]]<-weightstmp[kchoose[2]]*weightscomb*R1/weightsdash
#weights[kchoose[2]]<-weightstmp[kchoose[1]]*weightscomb*R2/weightsdash
#v[kchoose[2]]<-weights[kchoose[2]]/prod(1-v[1:(kchoose[2]-1)])
#v[kchoose[1]]<-weights[kchoose[1]]/prod(1-v[1:(kchoose[1]-1)])
#weights<-weightsupdate(Kstar,v)
#}
}
output<-list(z=z,nk=nk,v=v,weights=weights,mu=mu,Sigma=Sigma)
return(output)
}

updatealpha<-function(n,alpha,Kstar,a,b){
a<-1;b<-1;
for (t in 1:50){
phi1<-rbinom(1,1,n/(n+alpha))
phi<-rbeta(1,alpha+1,n)
alpha<-rgamma(1,shape=a+Kstar-phi1,rate=(b-log(phi)))
}
return(alpha)
}