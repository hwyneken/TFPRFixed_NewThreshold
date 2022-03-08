rm(list = ls())
require(ggplot2)
require(ggthemes)
require(scales)
require(tidyr)


Perc <-  0.5 #seq (0, 0.75, 0.25)
nRep <- 2000 # needs to be divisible by 4 
nTrial <- 5
n <- 100 #c(20, 50, 100, 400)
ClinSigLevels <- c(0,0.15,0.3,0.45) 
ClinSigLevels_Char <- paste0("ClinSig = ",ClinSigLevels)

### read in simulation data
EffectSizes <- readRDS(file = "sim_data/EffectPoint0.RDS")
EffectSizes <- as.numeric(EffectSizes)
sevScoresArray <- readRDS(file = "sim_data/sevTop2Array.RDS")
bfArray <- readRDS(file = "sim_data/BFPoint0t5.RDS")


### reorganize data

# re-use this function for the first element
getFirst <- function(x) {
  res <- head(x,1)
}

sevScores <- apply(sevScoresArray,c(3,4),
                   getFirst)
NSim <- dim(sevScores)[1]
colnames(sevScores) <- ClinSigLevels_Char
bfScores <- apply(bfArray,c(3,4),
                  getFirst)
logBFScores <- log(bfScores)
colnames(logBFScores) <- ClinSigLevels_Char


plotDF1 <- data.frame(EffectSizeCat = c(ifelse(EffectSizes > 0,"ClinSig","Below ClinSig"),
                                        ifelse(EffectSizes > 0.15,"ClinSig","Below ClinSig"),
                                        ifelse(EffectSizes > 0.3,"ClinSig","Below ClinSig"),
                                        ifelse(EffectSizes > 0.45,"ClinSig","Below ClinSig")),
                      EffectSize = rep(EffectSizes,times=4),
                      SeverityScore = c(sevScores[,1],sevScores[,2],
                                        sevScores[,3],sevScores[,4]),
                      BFScore = c(logBFScores[,1],logBFScores[,2],
                                  logBFScores[,3],logBFScores[,4]),
                      ClinSigThresh = rep(ClinSigLevels_Char,each=NSim))

plot1 <- ggplot(data=plotDF1,aes(x=SeverityScore,y=EffectSize))
plot1 <- plot1 + geom_point() 
plot1 <- plot1 + facet_wrap(ClinSigThresh ~ EffectSizeCat,nrow=4,
                            scales="free_x")
plot1 <- plot1 + theme_tufte()
plot1 <- plot1 + labs(x = "SEV Score",y="Standardized Effect Size",
                      title = "Changes in Severity Score with Effect Size\nAnd Clinical Relevance Threshold")
plot1 <- plot1 + theme(strip.background = element_rect(fill="gray"),
                       panel.background = element_rect(fill="#f8f8f8"),
                       strip.text = element_text(size=14),
                       axis.title = element_text(size=14),
                       axis.text = element_text(size=12),
                       title = element_text(size=16))
ggsave(plot1,file="images/SevVsEffectSize.png",units="in",
       width=12,height=8)


### Plot 2 - log pseudo BF vs Effect Size
plot2 <- ggplot(data=plotDF1,aes(x=BFScore,y=EffectSize))
plot2 <- plot2 + geom_point() 
plot2 <- plot2 + facet_wrap(ClinSigThresh ~ EffectSizeCat,nrow=4,
                            scales="free_x")
plot2 <- plot2 + theme_tufte()
plot2 <- plot2 + labs(x = "log(Pseudo-Bayes Factor) Score",y="Standardized Effect Size",
                      title = "Changes in Pseudo Bayes Factor with Effect Size\nAnd Clinical Relevance Threshold")
plot2 <- plot2 + theme(strip.background = element_rect(fill="gray"),
                       panel.background = element_rect(fill="#f8f8f8"),
                       strip.text = element_text(size=14),
                       axis.title = element_text(size=14),
                       axis.text = element_text(size=12),
                       title = element_text(size=16))
ggsave(plot2,file="images/BFVsEffectSize.png",units="in",
       width=12,height=8)

### Step 3 - create overall TP and FP rate graph

pvalueTP <- readRDS(file = "sim_data/pvalueTP.RDS")
pvalueFP <- readRDS(file = "sim_data/pvalueFP.RDS")
sevTop2TP <- readRDS(file=  "sim_data/sevTop2Hit.RDS")
sevTop2FP <- readRDS(file = "sim_data/sevTop2FA.RDS")
bfTP <- readRDS(file = "sim_data/bfTP.RDS")
bfFP <- readRDS(file = "sim_data/bfFP.RDS")

## collapse into numeric vector
pvalueTP <- apply(pvalueTP,3,getFirst)
pvalueFP <- apply(pvalueFP,3,getFirst)
sevTop2TP <- apply(sevTop2TP,3,getFirst)
sevTop2FP <- apply(sevTop2FP,3,getFirst)
bfTP <- apply(bfTP,3,getFirst)
bfFP <- apply(bfFP,3,getFirst)

# pull plot data together
tpVec <- c(pvalueTP,sevTop2TP,bfTP)
fpVec <- c(pvalueFP,sevTop2FP,bfFP)
csVec <- factor(rep(ClinSigLevels_Char,times=3),
                levels=ClinSigLevels_Char)
MethodVec <- factor(rep(c("P-Value Thresh",
                          "Severity Top-2",
                          "Pseudo BF > 10"),
                        each = length(pvalueTP)),
                    levels = c("P-Value Thresh",
                               "Severity Top-2",
                               "Pseudo BF > 10"))
plotDF3 <- data.frame(FP = fpVec,
                      TP = tpVec,
                      CS_Level = csVec,
                      Method = MethodVec)

plot3 <- ggplot(data=plotDF3,aes(x=FP,y=TP,color=Method,shape=Method))
plot3 <- plot3 + geom_point(size=6)
plot3 <- plot3 + facet_wrap(~ CS_Level,nrow=2,
                            scales="free_y")
plot3 <- plot3 + theme_tufte()
plot3 <- plot3 + labs(x = "False Positive Rate",
                      y = "True Positive Rate",
                      title = "False vs True Positive Clinical Trial Outcomes\nTPR = 50%, sample size/arm = 100",
                      caption = "Comparable to the Top Right of Figure 1 of Ravenzwaaij and Ioannidis 2019",
                      color="Method",shape="Method")
plot3 <- plot3 + theme(strip.background = element_rect(fill="gray"),
                       panel.background = element_rect(fill="#f8f8f8"),
                       strip.text = element_text(size=14),
                       axis.title = element_text(size=14),
                       axis.text = element_text(size=12),
                       title = element_text(size=16),
                       legend.title = element_text(size=14),
                       legend.text = element_text(size=12))
ggsave(plot3,file = "images/FPvsTP_Comp.png",units="in",
       width=12,height=8)

