### The code that makes this data is in "rand_power_indices.R". The simulations were run from the file "run_power_rand.R".
library(tidyverse)
library(gridExtra)
library(patchwork)

#naivezcol <- "#7CAE00"
#sampsplitcol <- "#00BFC4"
treevalcol <- "#F8766D"
sigmaSSEcol <- "navy"
sigmaConscol <- "darkgray"
#treevalESTcol <- "red"
#ctreecol <- "black"


setwd("~/Dropbox/Tree Values Paper/Code Reviewer Response/unknown variance/")
res<- read.csv("Full_Rand_Comps_4-8-2022.csv", header=FALSE,sep=" ")
names(res)<- c("cp", "seed", "depth", "beta", "method", "pval", 
               "pval_cons", "pval_SSE", "pval_med", "pval_mad", "truesig","n1", "n2", "samplesig", "bestRand", "whichSplit","XORlev", "best", "NA")



ggplot(data=res %>% filter(!is.na(depth), !is.na(XORlev)))+
  geom_smooth(aes(x=abs(truesig), y=as.numeric(pval < 0.05), col="known"), lwd=2, se=F, method="glm", method.args=list(family="binomial"))+
  geom_smooth(aes(x=abs(truesig), y=as.numeric(pval_cons < 0.05), col="cons"), se=F, method="glm", method.args=list(family="binomial"))+
  geom_smooth(aes(x=abs(truesig), y=as.numeric(pval_SSE < 0.05), col="SSE"), se=F, method="glm", method.args=list(family="binomial"))+
  theme_bw()+facet_grid(rows=vars(XORlev), col=vars(depth))



#res <- res %>% filter(!is.na(bestRand))
# res2 <- res %>% group_by(cp, method, seed,XORlev,beta, depth, pval, 
#                          pval_cons, pval_SSE, pval_med, pval_mad, truesig, n1, n2, samplesig) %>% 
#   summarize(bestRand = max(bestRand),whichSplit = which.max(bestRand))
res2 <- res 
res2$whichSplit <- as.numeric(res2$whichSplit)
res2 <- res2 %>% mutate(whichSplitLev = ifelse(whichSplit >=3,3,whichSplit))
res2 <- res2 %>% mutate(truesplitlev = ifelse(whichSplit >=3, 3, whichSplit))
res2 <- res2 %>%
  mutate(
    correct = as.numeric(bestRand > 0.75),
    correctreject = as.numeric(bestRand > 0.75 & pval < 0.05),
    correctreject_cons = as.numeric(bestRand > 0.75 & pval_cons < 0.05),
    correctreject_SSE = as.numeric(bestRand > 0.75 & pval_SSE < 0.05),
    correctreject_med = as.numeric(bestRand > 0.75 & pval_med < 0.05)
  )

#### HOW MANY TRUE TREES WERE MADE FOR EACH BETA. It's just 500, which I set. Didn't want it hardcoded though.
ntree <-  NROW((res %>% filter(depth==1,method=="condition",XORlev==2,beta==10)))


groupRes <- res2 %>% group_by(method, XORlev, truesplitlev, beta) %>% summarize(
  numcorrect = sum(correct),
  numcorrectreject = sum(correctreject),
  numcorrectreject_cons = sum(correctreject_cons),
  numcorrectreject_SSE = sum(correctreject_SSE),
  numcorrectreject_med = sum(correctreject_med)
) %>%
  mutate(propCorrect = ifelse(truesplitlev < 3, numcorrect/ntree, numcorrect/(2*ntree)),
         propCorrectReject= ifelse(truesplitlev < 3, numcorrectreject/ntree, numcorrectreject/(2*ntree)),
                                   propCorrectReject= ifelse(truesplitlev < 3, numcorrectreject/ntree, numcorrectreject/(2*ntree)),
                                   propCorrectReject_cons = ifelse(truesplitlev < 3, numcorrectreject_cons/ntree, numcorrectreject_cons/(2*ntree)),
                                   propCorrectReject_SSE = ifelse(truesplitlev < 3, numcorrectreject_SSE/ntree, numcorrectreject_SSE/(2*ntree)),
                                   propCorrectReject_med = ifelse(truesplitlev < 3, numcorrectreject_med/ntree, numcorrectreject_med/(2*ntree))

  ) %>%
  mutate(
    beta=as.numeric(beta),
    propCorrect=as.numeric(propCorrect),
    propCorrectReject=as.numeric(propCorrectReject),
    propCorrectReject_cons = as.numeric(propCorrectReject_cons)
  ) %>%
  mutate(truesplitlev2 = paste("Level ", truesplitlev))


groupRes$XORlev2 <- paste0("a=", groupRes$XORlev)

##### Estimated vs Not estimated. POwer ONly
groupRes$level2 <- paste("Level", groupRes$truesplitlev)
p3 <- ggplot(data=groupRes %>% filter(XORlev==1, !is.na(truesplitlev), !is.na(XORlev2))) +
  geom_smooth(aes(x=beta,y=propCorrect, col="zknown", lty="Proportion Detected"),
              se=F,method="glm",  method.args = list(family = "binomial"))+
  geom_smooth(aes(x=beta,y=propCorrectReject, col="zknown", lty="Proportion Detected and Rejected"),
              se=F,method="glm",  method.args = list(family = "binomial"))+
  geom_smooth(aes(x=beta,y=propCorrectReject_cons, col="cons",lty="Proportion Detected and Rejected"),
              se=F,method="glm",  method.args = list(family = "binomial"))+
  geom_smooth(aes(x=beta,y=propCorrectReject_SSE, col="SSE"),
              se=F,method="glm",  method.args = list(family = "binomial"), lty=3)+
  ylab("Proportion of True Splits Detected and Rejected")+xlab("b")+
  facet_grid(cols  = 
               vars(level2)) + theme_bw()+
  scale_color_manual(labels=c(expression(hat(sigma)[cons]), expression(hat(sigma)[SSE]),
                              expression(paste("True ", sigma))),
                       values=c(sigmaConscol, sigmaSSEcol, treevalcol))+
  theme(legend.text.align = 0)+
  labs(col="", linetype="")+
  guides(linetype="none")
p3
ggsave("~/Dropbox/Tree Values Paper/JMLR Reviewer Responses/power.png",
       width=8, height=4)
ggsave("~/Dropbox/Tree Values Paper/JMLR_submit_v3/Figures/power.png",
       width=8, height=4)


