#2_two_component_normal.R
#generates figures to show different solutions for a two component Normal mixture
library(tidyverse)
library(ggpubr)

set.seed(12345)
#first row of plots: keep mu's fixed but vary eta
eta_1 = c(0.1,0.5,0.8)

mu_1 = 5
mu_2 = 2
sd_1 = sd_2 = 1

#simulate data
n = 10000
M = length(eta_1)


dat = vector('list',M)
for (t in 1:M){
  #sample latent indicate
  z = sample(1:2,size=n,prob=c(eta_1[t],1-eta_1[t]),replace=T)
  y.mix = (z==1)*rnorm(n,mu_1,sd_1) + (z==2)*rnorm(n,mu_2,sd_2)
  dat[[t]] = tibble(z=factor(z),y=y.mix)
  
}
#bind_rows before plotting
dat = bind_rows(dat,.id='scenario')

#plot as smoothed densities
#scenario 1
temp_dat = filter(dat,scenario==1)
g1 = ggplot(temp_dat, aes(y)) + geom_density(size=2) +
  theme_minimal()+
  expand_limits(x=c(-10,10),y=c(0,0.4))+
  scale_x_continuous('y',breaks = seq(-10,10,2))+
  scale_y_continuous('Density')+theme(axis.text = element_text(size=12))

#scenario 2
temp_dat = filter(dat,scenario==2)
g2 = ggplot(temp_dat, aes(y)) + geom_density(size=2) +
  theme_minimal()+
  expand_limits(x=c(-10,10),y=c(0,0.4))+
  scale_x_continuous('y',breaks = seq(-10,10,2))+
  scale_y_continuous('Density')+theme(axis.text = element_text(size=12))


#scenario 3
temp_dat = filter(dat,scenario==3)
g3 = ggplot(temp_dat, aes(y)) + geom_density(size=2) +
  theme_minimal()+
  expand_limits(x=c(-10,10),y=c(0,0.4))+
  scale_x_continuous('y',breaks = seq(-10,10,2))+
  scale_y_continuous('Density')+theme(axis.text = element_text(size=12))

#g4,g5,g6 here
eta_1 = 0.2
mu_2 = c(-2,0,2)
dat = vector('list',M)
z = sample(1:2,size=n,prob=c(eta_1,1-eta_1),replace=T)
for (t in 1:M){
  y.mix = (z==1)*rnorm(n,mu_1,sd_1) + (z==2)*rnorm(n,mu_2[t],sd_2)
  dat[[t]] = tibble(z=factor(z),y=y.mix)
  
}
#bind_rows before plotting
dat = bind_rows(dat,.id='scenario')

#plot as smoothed densities
#scenario 1
temp_dat = filter(dat,scenario==1)
g4 = ggplot(temp_dat, aes(y)) + geom_density(size=2) +
  theme_minimal()+
  expand_limits(x=c(-10,10),y=c(0,0.4))+
  scale_x_continuous('y',breaks = seq(-10,10,2))+
  scale_y_continuous('Density')+theme(axis.text = element_text(size=12))

#scenario 2
temp_dat = filter(dat,scenario==2)
g5 = ggplot(temp_dat, aes(y)) + geom_density(size=2) +
  theme_minimal()+
  expand_limits(x=c(-10,10),y=c(0,0.4))+
  scale_x_continuous('y',breaks = seq(-10,10,2))+
  scale_y_continuous('Density')+theme(axis.text = element_text(size=12))


#scenario 3
temp_dat = filter(dat,scenario==3)
g6 = ggplot(temp_dat, aes(y)) + geom_density(size=2) +
  theme_minimal()+
  expand_limits(x=c(-10,10),y=c(0,0.4))+
  scale_x_continuous('y',breaks = seq(-10,10,2))+
  scale_y_continuous('Density')+theme(axis.text = element_text(size=12))



ggarrange(g1,g2,g3,g4,g5,g6,nrow=2,ncol=3,align='hv')
ggsave("figures/example.2componentnormal.png", width = 30, height = 10, units = "cm")
invisible(dev.off())
