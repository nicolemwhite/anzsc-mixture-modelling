#profile_regression_example.R
library(PReMiuM)


#This script shows examples of profile regression applied to simulated datasets in the PReMiuM package

# example for Poisson outcome and Discrete covariates
inputs <- generateSampleDataFile(clusSummaryPoissonDiscrete())

runInfoObj<-profRegr(yModel=inputs$yModel,
                     xModel=inputs$xModel, nSweeps=10, nClusInit=20,
                     nBurn=20, data=inputs$inputData, output="output",
                     covNames = inputs$covNames, outcomeT = inputs$outcomeT,
                     fixedEffectsNames = inputs$fixedEffectNames)

dissimObj<-calcDissimilarityMatrix(runInfoObj)
clusObj<-calcOptimalClustering(dissimObj)
riskProfileObj<-calcAvgRiskAndProfile(clusObj)
clusterOrderObj<-plotRiskProfile(riskProfileObj,"figures/profile.regression.poisson.dicreate.summary.png")
