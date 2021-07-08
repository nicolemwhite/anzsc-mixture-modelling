#clustering_inference.R
## This scripts covers approaches to deal with label switching and estimating likely clusterings using MCMC output
## R dependencies: mcclust, reshape2::melt, tidyverse
require(mcclust)
require(reshape2)
require(tidyverse)


#load results from fmm_mvn.R as an example
load('mcmc output/fmm_mvn_MCMCoutput.rda')


#label switching approaches applied to sample MCMC output

## apply Stephens' relabelling algorithm. Uses the relabel() function
## outputs from this function include 
### unswitched labels over MCMC output (z.mcmc.post)
### posterior probabilities of membership based on relablled output (zprob.post). Really useful for looking into classification uncertainty
### most likely classification (z.post)

z_unswitched = relabel(t(MCMC.traces$z))
z.mcmc.post = z_unswitched$cls 
zprob.post = z_unswitched$P
z.post = z_unswitched$cl

### plot the posterior probabilities of membership to look at the uncertainty
zprob.post_long = reshape2::melt(zprob.post,varnames = c('obs','cluster'),value.name='Posterior probability')

#plot
ggplot(zprob.post_long,aes(x=cluster,y=obs,fill=`Posterior probability`))+
  geom_raster()+scale_x_continuous('Cluster')+
  scale_y_continuous('Observation index')

ggsave("figures/fmm_mvn_posteriorprobs.png", width = 25, height = 20, units = "cm")
invisible(dev.off())

## similarity matrix approaches (focusses on similarities in cluster assignment)

## apply similarity matrix to visualise clustering
z.sim = comp.psm(t(MCMC.traces$z))

#reorder for heatmap (easier to visualise)
tmp<-hclust(dist(z.sim))
z.sim_reorder<-z.sim[tmp$order,tmp$order]
z.sim_long <- reshape2::melt(z.sim_reorder,varnames = c('z1','z2'),value.name='Posterior probability')

#plot
ggplot(z.sim_long,aes(x=z1,y=z2,fill=`Posterior probability`))+
  geom_raster()+scale_x_continuous('Observation index',breaks=seq(0,350,25))+
  scale_y_continuous('Observation index',breaks=seq(0,350,25))

ggsave("figures/fmm_mvn_similarity.png", width = 25, height = 20, units = "cm")
invisible(dev.off())
