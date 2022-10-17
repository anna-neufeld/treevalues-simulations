library(tidyverse)
library(gridExtra)

naivezcol <- "#7CAE00"
sampsplitcol <- "#00BFC4"
treevalcol <- "#F8766D"

setwd("~/Dropbox/Tree Values Paper/Code Reviewer Response/non normal/")
nullResPoisson <- read.csv("poisson.csv", sep=" ", header=FALSE)
nullResGamma <- read.csv("gamma.csv", sep=" ", header=FALSE)
nullResBernoulli <- read.csv("bernoulli.csv", sep=" ", header=FALSE)
nullResCauchy <- read.csv("cauchy.csv", sep=" ", header=FALSE)
nullResMix <- read.csv("mixture.csv", sep=" ", header=FALSE)
nullResBernSmall <- read.csv("bernoulli01.csv", sep=" ", header=FALSE)
dat <- rbind(
  cbind(nullResPoisson, "dist"="Poisson(10)"),
  cbind(nullResGamma, "dist"= "Gamma(1,10)"),
  cbind(nullResBernoulli, "dist"="Bernoulli(0.5)"),
  cbind(nullResCauchy, "dist" = "Cauchy(0,1)"),
  cbind(nullResMix, "dist" = "Mixture"),
  cbind(nullResBernSmall, "dist" = "Bernoulli(0.1)")
)
names(dat) <- c("cp", "seed", "method", "type", "depth", "pval", "n1", "n2", "sampsig", "NA", "dist")

dat$method[dat$method=="splitting"] = "Sample Splitting"
dat$method[dat$method=="NaiveZ"] = "Naive Z-test"
dat$method[dat$method=="Tree-Values"] = "Selective Z-test"
dat$method <- ordered(dat$method, levels=c("Naive Z-test", "Sample Splitting", "Selective Z-test"))
dat$dist <- ordered(dat$dist, levels=c("Poisson(10)", "Bernoulli(0.5)", "Gamma(1,10)","Cauchy(0,1)",
                                       "Miture", "Bernoulli(0.1)"))
ggplot(data=dat %>% filter(dist != "Mixture",
                           dist != "Cauchy(0,1)"), aes(sample=pval, col=method))+geom_qq(distribution="qunif")+
  facet_wrap(vars(dist))+coord_fixed()+theme_bw()+
  geom_abline(slope=1,intercept=0)+
  labs(col="Method")+
  scale_color_manual(values = c(naivezcol, sampsplitcol, treevalcol))+
  xlab("Unif(0,1) Quantiles")+ylab("Empirical Quantiles")

ggsave("~/Dropbox/Tree Values Paper/JMLR Reviewer Responses/non_normal.png")
ggsave("~/Dropbox/Tree Values Paper/JMLR_submit_v3/Figures/non_normal.png")
