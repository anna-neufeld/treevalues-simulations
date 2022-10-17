library(tidyverse)
library(gridExtra)
library(patchwork)

setwd("~/Dropbox/Tree Values Paper/Code Reviewer Response/unknown variance")
res <- read.csv("non_null_var_est_4-5-2022.csv", header=FALSE, sep=" ")
names(res) <- c("beta", "XORLev", "cp", "seed", "method", "type", "depth", 
                "pval", "pval_cons", "pval_SSE", "pval_MAD", "pval_MED",
                "lower", "upper",  "lower_cons", "upper_cons",  "lower_SSE", "upper_SSE",  "lower_MAD", "upper_MAD",  "lower_MED", "upper_MED", 
                "n1", "n2", "sampsig","truesig", "truesigtest", "sizeS", "sizeNum")
res <- res %>% mutate(correctCI = lower < truesig & upper > truesig, 
                      correctCI_cons= lower_cons < truesig & upper_cons > truesig,
                      correctCI_SSE= lower_SSE < truesig & upper_SSE > truesig,
                      correctCI_mad= lower_MAD < truesig & upper_MAD > truesig,
                      correctCI_med= lower_MED < truesig & upper_MED > truesig,
                      length = upper - lower,
                      length_cons = upper_cons-lower_cons,
                      length_SSE = upper_SSE-lower_SSE)


###  TABLE 2 !!!
res %>% group_by(method, type, depth) %>% summarize("Known" = round(mean(correctCI, na.rm=TRUE),3),
                                                    "Cons" = round(mean(correctCI_cons, na.rm=TRUE),3),
                                                    "SSE" = round(mean(correctCI_SSE, na.rm=TRUE),3))

ggplot(data=res %>% filter(truesig==0))+
  geom_qq(aes(sample=pval_cons, col="Conservative"), distribution="qunif")+
  geom_qq(aes(sample=pval_SSE, col="SSE"), distribution="qunif")+
  geom_qq(aes(sample=pval, col="zKnown"), distribution="qunif")+
  geom_abline(slope=1,intercept=0)+theme_bw()+
  xlab("Unif(0,1) Quantiles")+ylab("Sample Quantiles")+
  facet_grid(rows=vars(type), cols=vars(depth))+
  scale_color_discrete(labels=c(expression(hat(sigma)[cons]), expression(hat(sigma)[SSE]),
                              expression(sigma)))+
  labs(col="Variance Estimate")+coord_fixed()
ggsave("~/Dropbox/Tree Values Paper/JMLR Reviewer Responses/null_fig.png",
       width=8, height=4)
