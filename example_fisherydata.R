library(bayesmix)
library(tidyverse)

data("fish", package = "bayesmix")

ggplot(fish,aes(fish))+geom_histogram(bins=40,fill='white',colour='black')+
theme_minimal()+scale_x_continuous('length in inches')
ggsave("figures/fishery.data.png", width = 25, height = 20, units = "cm")
invisible(dev.off())


model <- BMMmodel(fish, k = 4, priors = list(kind = "independence",
                                             parameter = "priorsFish", hierarchical = "tau"),
                  initialValues = list(S0 = 2))
model
