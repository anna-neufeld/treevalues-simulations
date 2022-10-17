setwd("~/Dropbox/Tree Values Paper/Code Reviewer Response/more permutations/")
library(tidyverse)
filename <- "permutations_3-32-2022.csv"
res <- read.csv(filename, header=FALSE, sep=" ")
names(res) <- c("beta","XORlev", "seed", "depth", "samp_sig", "true_sig",
                  "lengths01", "lengths1", "pvals1", "CIs1_low","CIs1_high",
                  "lengths02", "lengths2", "pvals2","CIs2_low","CIs2_high",
                  "lengths03", "lengths3", "pvals3","CIs3_low","CIs3_high",
                  "lengths04", "lengths4", "pvals4","CIs4_low","CIs4_high",
                  "lengths05", "lengths5", "pvals5","CIs5_low","CIs5_high",
                  "lengths06", "lengths6", "pvals6","CIs6_low","CIs6_high")

res <- res %>% mutate(width1 = CIs1_high-CIs1_low,
                      width2 = CIs2_high-CIs2_low,
                      width3 = CIs3_high-CIs3_low,
                      width4 = CIs4_high-CIs4_low,
                      width5 = CIs5_high-CIs5_low,
                      width6 = CIs6_high-CIs6_low)

res <- res %>% filter(XORlev==1)

#### ITS NOT FOR SIBS so I CANT DO THAT DETECTION DEFINITION
#### I THINK this conditional power plot makes the most sense
p1 <- ggplot(data=res)+geom_smooth( aes(x=true_sig, y=as.numeric(pvals1 < 0.05), col="Identity permutation only"), se=F, method="glm",  method.args = list(family = "binomial"))+
  #geom_smooth(aes(x=true_sig, y=as.numeric(pvals2 < 0.05), col="2"), se=F, method="glm",  method.args = list(family = "binomial"))+
  #geom_smooth( aes(x=true_sig, y=as.numeric(pvals3 < 0.05), col="3"), se=F, method="glm",  method.args = list(family = "binomial"))+
  #geom_smooth( aes(x=true_sig, y=as.numeric(pvals4 < 0.05), col="4"), se=F, method="glm",  method.args = list(family = "binomial"))+
  #geom_smooth( aes(x=true_sig, y=as.numeric(pvals5 < 0.05), col="5"), se=F, method="glm",  method.args = list(family = "binomial"))+
  geom_smooth( aes(x=true_sig, y=as.numeric(pvals6 < 0.05), col= "Full conditioning set"), se=F, method="glm",  method.args = list(family = "binomial"))+
  xlim(0,20)+
  ylab("Estimated Power")+xlab("True Difference in Means")+labs(col="Method")+
  theme_bw()#+facet_grid(vars(XORlev))


consRes <- res %>% group_by(beta, XORlev) %>% summarize(width1=median(width1),
                                                        width2 = median(width2),
                                                        width3 = median(width3),
                                                        width4 = median(width4),
                                                        width5 = median(width5),
                                                        width6 = median(width6)) 
p2 <- ggplot(data=consRes)+geom_smooth( aes(x=beta, y=width1, col="Identity permutation only"), se=F)+
  #geom_smooth(aes(x=beta, y=width2, col="2"), se=F)+
  #geom_smooth( aes(x=beta, y=width3, col="3"), se=F)+
  #geom_smooth( aes(x=beta, y=width4, col="4"), se=F)+
  #geom_smooth( aes(x=beta, y=width5, col="5"), se=F)+
  geom_smooth( aes(x=beta, y=width6, col="Full conditioning set"), se=F)+
  xlab("b")+ylab("Median CI Width") + theme_bw()+labs(col="Method") #+facet_grid(vars(XORlev))

p3 <- p1 + xlim(7.5,15)
library(patchwork)
p1+p3+p2+plot_layout(guides="collect")
ggsave("~/Dropbox/Tree Values Paper/JMLR Reviewer Responses/permutation_res.png", width=9, height=3)
ggsave("~/Dropbox/Tree Values Paper/JMLR_submit_v3/Figures/permutation_res.png", width=9, height=3)


