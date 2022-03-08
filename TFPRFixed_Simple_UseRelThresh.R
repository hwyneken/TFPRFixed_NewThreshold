rm (list = ls())
library (BayesFactor)
library (pCalibrate)  # MinBF01, or max evidence against H0
set.seed(7788)

## CASE References
# Original Paper: https://bmcmedresmethodol.biomedcentral.com/articles/10.1186/s12874-019-0865-y
# OSF link to original code: https://osf.io/t2dev/
# Order statistics for the uniform: https://en.wikipedia.org/wiki/Order_statistic#Order_statistics_sampled_from_a_uniform_distribution
# much of the code/variable name conventions are copied directly from the original
# COS ToS https://github.com/CenterForOpenScience/cos.io/blob/master/TERMS_OF_USE.md

#### Simulate Clinical Trial Outcomes
############################################

### Input
# We don't need to load in any data files
# The main input step is defining the simulation parameters
# Params: 
#   Perc: The probability that the experimental group actually has a different mean
#     e.g. at 0.5, for 50% of the simulations, the experimental group delta = 0
#   nRep: the number of indep simulations (2000)
#   n: The sample size per arm (100)
#   ClinSigLevels: the standardized exp vs control difference deemed to be meaningful by the FDA
#   
#   We will also allocate arrays to hold the results
#   EffectPoint0: holds the fixed-effect sizes for each simulation
#   pPoint0: holds the one-sided (exp > ctrl) p-values each trial within each simulation rep
#   sevTop2Array: holds the severity scores (def below) for each sim rep and clinsig level
#   BFPoint0t5: holds the pseduo-Bayes scores for each sim rep and clinsig level

### Output
# 
# the script ends by storing the above-mentioned arrays, plus summaries of
#   true and false positives according to their decision rules


### Details
#
# Algorithm: 
#   1) cycle through all choices of Perc, n and sim rep
#     generate effect sizes and data
#   2) start an inner loop for choices of clinical significance
#     for example, calculate p-values using a null difference = 0, 0.15, etc...
#     within this inner loop, calculate p-values, severity scores and pseudo bayes factors
#   3) calculate true/false positive rates
#     use the following decision rules: 
#       p-values: claim discovery if > 1/5 trials has a p-value < 0.025
#       severity: claim discovery if severity >= 0.95
#         using a threshold of pbinom(1,5,0.025) gave the exact same results as the p-values
#       pseudo-BF: claim discovery if BF >= 10
# Definitions
#   Severity: taken from Mayo 2018 
#     somewhat arbitrarily choose 
#     SEV = 1 - P(p(2) < p(2)_Obs ; H0)
#       where p(2) is the 2nd order statistic of the five p-values - has Beta(2,4) dist


Perc <-  0.5 #seq (0, 0.75, 0.25)
nRep <- 2000 # needs to be divisible by 4 
nTrial <- 5
n <- 100 #c(20, 50, 100, 400)
EffectPoint0 <- matrix (, length (Perc), nRep)
ClinSigLevels <- c(0,0.15,0.3,0.45) 

pPoint0 <- array (, dim = c(length (Perc), length (n), nRep, nTrial,length(ClinSigLevels)))
sevTop2Array <- array (, dim = c(length (Perc), length (n), nRep, length(ClinSigLevels)))
BFPoint0t5 <- array (, dim = c(length (Perc), length (n), nRep, length(ClinSigLevels)))

### Step 1

for (hh in 1:length(Perc)) { # HE
  EffectPoint0[hh,] <- c(rep (0, nRep*Perc[hh]), rnorm (nRep*(1-Perc[hh]), 0.4, 0.13))
  for (h in 1:length (n)) {
    for (i in 1:nRep) {
      cat ("\n ", paste (hh, "of", length (Perc), ",", h, "of", length (n), ",", i, "of", nRep, sep = " "))
      Plac <- matrix (rnorm (nTrial*n[h], 0, 1), nTrial, n[h])
      Treat <- matrix (rnorm (nTrial*n[h], EffectPoint0[hh,i], 1), nTrial, n[h])

      
      for (m in 1:nTrial) {
        for (m2 in 1:(length(ClinSigLevels))) {
          pPoint0[hh,h,i,m,m2] <- t.test (Treat[m,], Plac[m,], 
                                         alternative = "greater",
                                         mu=ClinSigLevels[m2])$p.value
        }
      }
      
      ### Step 2
      for (m in 1:(length(ClinSigLevels))) {
        tempBFRes <- as.vector (ttestBF (x = as.vector (Treat[1:5,]) - ClinSigLevels[m], 
                                                   y = as.vector (Plac[1:5,]), 
                                                   nullInterval = c(0, Inf)))[1]
        BFPoint0t5[hh,h,i,m] <- as.numeric(tempBFRes)
        
        
        tempPValues <- pPoint0[hh,h,i,,m]
        temp2ndSmallestPValue <- sort(tempPValues,
                                      decreasing=FALSE)[2]
        sevTop2Array[hh,h,i,m] <- 1 - pbeta(temp2ndSmallestPValue,
                                            shape1=2,shape2=4)
      }
      
    }
  }
}


#### Step 3
Alp <- 0.025
pHitPoint0t5 <- array (, dim = c(length (Perc), length (n), length (ClinSigLevels)))
pFAPoint0t5 <- array (, dim = c(length (Perc), length (n), length (ClinSigLevels))) 

SevThresh <- 0.95 # pbinom(1,size=5,prob=0.025)
sevTop2Hit <- array (, dim = c(length (Perc), length (n), length (ClinSigLevels)))
sevTop2FA <- array (, dim = c(length (Perc), length (n), length (ClinSigLevels)))

BFThresh <- 10 
BFHitPoint0t5 <- array (, dim = c(length (Perc), length (n), length (ClinSigLevels)))
BFFAPoint0t5 <- array (, dim = c(length (Perc), length (n), length (ClinSigLevels)))
  
for (hh in 1:(length(Perc))) {
  for (h in 1:(length(n))) {
    for (i in 1:(length(ClinSigLevels))) {
      
      ######### P-Value Cutoff
      tempTTestRes <- pPoint0[hh,h,,1:5,i]
      tempNumSigRes <- apply(tempTTestRes,1,
                             function (x) sum (x<Alp))
      tempClaimedDiscovery <- tempNumSigRes > 1
      
      tempMeaningfulEffects <- EffectPoint0[hh,] > ClinSigLevels[i]
      tempIrrelevantEffects <- EffectPoint0[hh,] <= ClinSigLevels[i]
      
      tempTruePositives <- tempClaimedDiscovery & tempMeaningfulEffects
      tempFalsePositives <- tempClaimedDiscovery & tempIrrelevantEffects
      
      tempNumME <- sum(tempMeaningfulEffects)
      tempNumIE <- sum(tempIrrelevantEffects)
      
      pHitPoint0t5[hh,h,i] <- ifelse(tempNumME == 0,NA,
                                     sum(tempTruePositives)/tempNumME)
      pFAPoint0t5[hh,h,i] <- ifelse(tempNumIE == 0,NA,
                                    sum(tempFalsePositives)/tempNumIE)
      
      #### Severity Cutoff
      tempSevRes <- sevTop2Array[hh,h,,i]
      tempClaimedDiscovery <- tempSevRes > SevThresh
      
      tempMeaningfulEffects <- EffectPoint0[hh,] > ClinSigLevels[i]
      tempIrrelevantEffects <- EffectPoint0[hh,] <= ClinSigLevels[i]
      
      tempTruePositives <- tempClaimedDiscovery & tempMeaningfulEffects
      tempFalsePositives <- tempClaimedDiscovery & tempIrrelevantEffects
      
      tempNumME <- sum(tempMeaningfulEffects)
      tempNumIE <- sum(tempIrrelevantEffects)
      
      sevTop2Hit[hh,h,i] <- ifelse(tempNumME == 0,NA,
                                     sum(tempTruePositives)/tempNumME)
      sevTop2FA[hh,h,i] <- ifelse(tempNumIE == 0,NA,
                                    sum(tempFalsePositives)/tempNumIE)
      
      
      ### BF Cutoff
      tempBFRes <- BFPoint0t5[hh,h,,i]
      
      tempClaimedDiscovery <- tempBFRes > BFThresh
      
      tempMeaningfulEffects <- EffectPoint0[hh,] > ClinSigLevels[i]
      tempIrrelevantEffects <- EffectPoint0[hh,] <= ClinSigLevels[i]
      
      tempTruePositives <- tempClaimedDiscovery & tempMeaningfulEffects
      tempFalsePositives <- tempClaimedDiscovery & tempIrrelevantEffects
      
      tempNumME <- sum(tempMeaningfulEffects)
      tempNumIE <- sum(tempIrrelevantEffects)
      
      BFHitPoint0t5[hh,h,i] <- ifelse(tempNumME == 0,NA,
                                      sum(tempTruePositives)/tempNumME)
      BFFAPoint0t5[hh,h,i] <- ifelse(tempNumIE == 0,NA,
                                     sum(tempFalsePositives)/tempNumIE)
    }
  }
}

saveRDS(EffectPoint0,file="sim_data/EffectPoint0.RDS")
saveRDS(sevTop2Array,file="sim_data/sevTop2Array.RDS")
saveRDS(BFPoint0t5,file="sim_data/BFPoint0t5.RDS")

saveRDS(pHitPoint0t5,file="sim_data/pvalueTP.RDS")
saveRDS(pFAPoint0t5,file="sim_data/pvalueFP.RDS")
saveRDS(sevTop2Hit,file="sim_data/sevTop2Hit.RDS")
saveRDS(sevTop2FA,file="sim_data/sevTop2FA.RDS")
saveRDS(BFHitPoint0t5,file="sim_data/bfTP.RDS")
saveRDS(BFFAPoint0t5,file="sim_data/bfFP.RDS")
