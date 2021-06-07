### The code that makes this data is in "rand_power_indices.R". The simulations were run from the file "run_power_rand.R".
library(tidyverse)
library(gridExtra)

naivezcol <- "#7CAE00"
sampsplitcol <- "#00BFC4"
treevalcol <- "#F8766D"
ctreecol <- "black"

setwd("~/Dropbox/Tree Values Paper/Code : Other for Creating All Figures")
res<- read.csv("Full_Rand_Comps_1-18-2021-maxdepth3.csv", header=FALSE,sep=" ")
names(res)<- c("cp", "seed", "depth", "beta", "method", "pval", "truesig","n1", "n2", "samplesig", "bestRand", "whichSplit","XORlev", "best", "NA")


res2 <- res %>% group_by(cp, method, seed,XORlev,beta, depth, pval, truesig, n1, n2, samplesig) %>% summarize(bestRand = max(bestRand), whichSplit = whichSplit[which.max(bestRand)])
res2 <- res2 %>% mutate(whichSplitLev = ifelse(whichSplit >=3,3,whichSplit))
res2 <- res2 %>% mutate(truesplitlev = ifelse(whichSplit >=3, 3, whichSplit))
res2 <- res2 %>%
  mutate(
    correct = as.numeric(bestRand > 0.75),
    correctreject = as.numeric(bestRand > 0.75 & pval < 0.05)
  )

#### HOW MANY TRUE TREES WERE MADE FOR EACH BETA. It's just 500, which I set. Didn't want it hardcoded though.
ntree <-  NROW((res %>% filter(depth==1,method=="condition",XORlev==2,beta==10)))


groupRes <- res2 %>% group_by(method, XORlev, truesplitlev, beta) %>% summarize(
  numcorrect = sum(correct),
  numcorrectreject = sum(correctreject)
) %>%
  mutate(propCorrect = ifelse(truesplitlev < 3, numcorrect/ntree, numcorrect/(2*ntree)),
         propCorrectReject= ifelse(truesplitlev < 3, numcorrectreject/ntree, numcorrectreject/(2*ntree))
  ) %>%
  mutate(
    beta=as.numeric(beta),
    propCorrect=as.numeric(propCorrect),
    propCorrectReject=as.numeric(propCorrectReject)
  ) %>%
  mutate(truesplitlev = paste("Level ", truesplitlev))


###  MAIN PLOT
groupRes$method <- factor(groupRes$method, ordered=TRUE, levels=c("splitting","ctree","condition"))
groupRes$XORlev2 <- paste0("a=", groupRes$XORlev)
p2 <- ggplot(data=groupRes %>% filter(beta > 0)
             , aes(x=beta, y=propCorrect, col=interaction(method, "detect"), lty=interaction(method, "detect"))) + geom_smooth(se=F,method="glm",  method.args = list(family = "binomial"))+
  geom_smooth(data=groupRes %>% filter(beta>0), aes(x=beta,y=propCorrectReject, col=interaction(method,"reject"), lty=interaction(method,"reject")),
              se=F,method="glm",  method.args = list(family = "binomial"))+
  ylab("Proportion of True Splits Detected or Rejected")+xlab("b")+
  scale_color_manual(name = "Method", values = c(treevalcol, treevalcol,ctreecol, ctreecol,sampsplitcol, sampsplitcol),
                     labels=c("CART with selective Z-tests, detection",
                              "CART with selective Z-tests, rejection",
                              "CTree, detection and rejection",
                              "",
                              "CART with sample splitting, detection",
                              "CART with sample splitting, rejection")) +
  scale_linetype_manual(name = "Method", values = c(1,3,1,0,1,3),
                        labels=c("CART with selective Z-tests, detection",
                                 "CART with selective Z-tests, rejection",
                                 "CTree, detection and rejection",
                                 "",
                                 "CART with sample splitting, detection",
                                 "CART with sample splitting, rejection"))+
  facet_grid(cols = vars(XORlev2), rows = vars(truesplitlev)) + theme_bw()+
  theme(legend.position=" ")

p2

#ggsave(filename="~/Dropbox/Tree Values Paper/newformat_v31/Figures/mainPower.png", p2, units="cm", width=16, height=11)
