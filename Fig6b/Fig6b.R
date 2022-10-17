library(tidyverse)
library(gridExtra)
library(patchwork)

naivezcol <- "#7CAE00"
sampsplitcol <- "#00BFC4"
treevalcol <- "#F8766D"
ctreecol <- "black"

setwd("~/Dropbox/Tree Values Paper/Code Reviewer Response/CI_RAND/")
res <- read.csv("non_null_var_est_10-10-2022.csv", header=FALSE, sep=" ")
names(res) <- c("beta", "XORLev", "cp", "seed", "method", "type", "depth", 
                "pval", 
                "lower", "upper", 
                "n1", "n2", "sampsig","truesig", "truesigtest", "sizeS", "sizeNum",
                "rand1", "rand2","rand3")
res <- res %>% mutate(correctCI = lower < truesig & upper > truesig, 
           length = upper - lower)
res$thisRand <- NA
res$thisRand[res$depth==1] <- res$rand1[res$depth==1]
res$thisRand[res$depth==2] <- res$rand2[res$depth==2] 
res$thisRand[res$depth==3] <- res$rand2[res$depth==3]


#### All levels, relate to your own levels RAND.
#### This is MEAN RAND. Very pessimistic

ggplot(data=res %>% filter(beta>0, type=="child"), aes(x=thisRand, y=length, group=as.factor(depth), col=method, lty=as.factor(depth)))+
  geom_quantile(quantiles=c(0.5),
                formula=y~splines::bs(x,df=3), lwd=1)+
  xlab("Adjusted Rand Index at Each Level")+ylab("Median Confidence Interval Width")+
  ylab("Median Confidence Interval Width")+xlab("Adjusted Rand Index (true vs. estimated tree)")+
  theme_bw()+
  scale_linetype_manual(values=c(1,2,3))+
  theme(axis.title=element_text(size=12), 
        title=element_text(size=12), axis.text=element_text(size=12))+
  guides(lty="none", col="none")
ggsave("~/Dropbox/Tree Values Paper/JMLR_resubmit_Oct_2022/Figures/rand_CI_median.png")


