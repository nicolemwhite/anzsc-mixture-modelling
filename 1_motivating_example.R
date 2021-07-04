#1_motviation_example.R
library(tidyverse)

m1 = 27
s1 = 2
m2 = 34
s2 = 4

n = 10000
eta = 0.6

z = sample(1:2,size=n,prob=c(eta,1-eta),replace=T)

y.mix = (z==1)*rnorm(n,m1,s1) + (z==2)*rnorm(n,m2,s2)

dat = tibble(z=factor(z),y=y.mix)


#overall distribution, no classification
ggplot(dat, aes(y)) + 
  geom_histogram(aes(y = stat(density)),binwidth=1,alpha = 0.3,colour='black') + theme_minimal()+
  scale_x_continuous('Body Mass Index (BMI)',breaks=seq(15,50,5))+
  scale_y_continuous('Density')+theme(axis.text = element_text(size=12))

ggsave("figures/motivation.null.png", width = 20, height = 20, units = "cm")
invisible(dev.off())

ggplot(dat, aes(y)) + 
  geom_histogram(aes(y = stat(density)),binwidth=1,alpha = 0.3,colour='black') + 
  stat_function(
    fun = dnorm, 
    args = list(mean = mean(dat$y), sd = sd(dat$y)), 
    lwd = 2, 
    col = 'red'
  )+
  theme_minimal()+
  scale_x_continuous('Body Mass Index (BMI)',breaks=seq(15,50,5))+
  scale_y_continuous('Density')+theme(axis.text = element_text(size=12))

ggsave("figures/motivation.null.singledist.png", width = 20, height = 20, units = "cm")
invisible(dev.off())


#known classification
ggplot(dat, aes(y)) + 
  geom_histogram(data = filter(dat,z=='1'), aes(y = stat(density)),fill = "#CC79A7", alpha = 0.3,binwidth=1,colour="#CC79A7") +
  stat_function(
    fun = dnorm, 
    args = list(mean = mean(filter(dat,z=='1')$y), sd = sd(filter(dat,z=='1')$y)), 
    lwd = 2, 
    col = "#CC79A7"
  )+ 
  geom_histogram(data = filter(dat,z=='2'), aes(y = stat(density)),fill = "#56B4E9", alpha = 0.3,binwidth=1,colour="#56B4E9") +
  stat_function(
    fun = dnorm, 
    args = list(mean = mean(filter(dat,z=='2')$y), sd = sd(filter(dat,z=='2')$y)), 
    lwd = 2, 
    col = "#56B4E9"
  )+   
  theme_minimal()+
  scale_x_continuous('Body Mass Index (BMI)',breaks=seq(15,50,5))+
  scale_y_continuous('Density')+theme(axis.text = element_text(size=12))
ggsave("figures/motivation.z.png", width = 20, height = 20, units = "cm")
invisible(dev.off())
