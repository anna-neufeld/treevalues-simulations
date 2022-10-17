library(tidyverse)
library(gridExtra)

naivezcol <- "#7CAE00"
sampsplitcol <- "#00BFC4"
treevalcol <- "#F8766D"
treevalESTcol <- "red"
ctreecol <- "black"

setwd("~/Dropbox/Tree Values Paper/Code Reviewer Response")
nullResNew <- read.csv("normal_var_est.csv", sep=" ", header=FALSE)
nullRes <- read.csv("null_res_1-18-2021-5000.csv",sep=" ", header=FALSE) 
names(nullRes) <- 
  names(nullResNew) <- c("cp", "seed", "method", "type", "depth", "pval", "n1", "n2", "sampsig", "NA")

quants_treevalue_est1 <- sort((nullResNew %>% filter(depth==1,method=="Tree-Values", type=="split"))$pval)
quants_treevalue_split1 <- sort((nullRes %>% filter(depth==1,method=="Tree-Values", type=="split"))$pval)
quants_pvalue_split1 <- sort((nullRes %>% filter(depth==1,method=="NaiveZ",type=="split"))$pval)
quants_sampsplit_split1 <- sort((nullRes %>% filter(depth==1,method=="splitting", type=="split"))$pval)

quants_treevalue_est2 <- sort((nullResNew %>% filter(depth==2,method=="Tree-Values", type=="split"))$pval)
quants_treevalue_split2 <- sort((nullRes %>% filter(depth==2,method=="Tree-Values", type=="split"))$pval)
quants_pvalue_split2 <- sort((nullRes %>% filter(depth==2,method=="NaiveZ",type=="split"))$pval)
quants_sampsplit_split2 <- sort((nullRes %>% filter(depth==2,method=="splitting", type=="split"))$pval)

quants_treevalue_est3 <- sort((nullResNew %>% filter(depth==3,method=="Tree-Values", type=="split"))$pval)
quants_treevalue_split3 <- sort((nullRes %>% filter(depth==3,method=="Tree-Values", type=="split"))$pval)
quants_pvalue_split3 <- sort((nullRes %>% filter(depth==3,method=="NaiveZ",type=="split"))$pval)
quants_sampsplit_split3 <- sort((nullRes %>% filter(depth==3,method=="splitting", type=="split"))$pval)

p1 <- ggplot(data=NULL) +
  geom_line(aes(x=qunif(seq(0,1,length.out=length(quants_pvalue_split1))), y=quants_pvalue_split1, col="Naive Z-test"), alpha=1, lwd=1.4) +
  geom_line(aes(x=qunif(seq(0,1,length.out=length(quants_sampsplit_split1))), y=quants_sampsplit_split1, col="Sample Splitting"),alpha=1, lwd=1.6)+
  geom_line(aes(x=qunif(seq(0,1,length.out=length(quants_treevalue_split1))), y=quants_treevalue_split1, col="Selective Z-test, Known Var"),alpha=1, lwd=1.1)+
  geom_line(aes(x=qunif(seq(0,1,length.out=length(quants_treevalue_est1))), y=quants_treevalue_est1, col="Selective Z-test, Estimated Var"),alpha=1, lwd=1.1)+
  xlab("") +
  ylab("Empirical Quantiles")+ggtitle("Level 1") +
  theme(plot.title = element_text(size = 10),axis.title=element_text(size = 7.5)) +
  coord_fixed()+theme_bw() + labs(col="Method")+guides(col=FALSE)+
  scale_color_manual(values = c(naivezcol, sampsplitcol, treevalcol, treevalESTcol))+
  geom_abline()

p2 <- ggplot(data=NULL) +
  geom_line(aes(x=qunif(seq(0,1,length.out=length(quants_pvalue_split2))), y=quants_pvalue_split2, col="Naive Z-test"), alpha=1, lwd=1.4) +
  geom_line(aes(x=qunif(seq(0,1,length.out=length(quants_sampsplit_split2))), y=quants_sampsplit_split2, col="Sample Splitting"),alpha=1, lwd=1.6)+
  geom_line(aes(x=qunif(seq(0,1,length.out=length(quants_treevalue_split2))), y=quants_treevalue_split2, col="Selective Z-test, Known Var"), alpha=1, lwd=1.1) +
  geom_line(aes(x=qunif(seq(0,1,length.out=length(quants_treevalue_est2))), y=quants_treevalue_est2, col="Selective Z-test, Estimated Var"),alpha=1, lwd=1.1)+
  xlab("U(0,1) Quantiles") +
  ylab("")+ggtitle("Level 2") +
  theme(plot.title = element_text(size = 10),axis.title=element_text(size = 7.5)) +
  coord_fixed()+theme_bw() + labs(col="Method")+guides(col=FALSE)+
  scale_color_manual(values = c(naivezcol, sampsplitcol, treevalcol, treevalESTcol))+
  geom_abline()

p3 <- ggplot(data=NULL)+
  geom_line(aes(x=qunif(seq(0,1,length.out=length(quants_pvalue_split3))), y=quants_pvalue_split3, col="Naive Z-test"), alpha=1, lwd=1.4) +
  geom_line(aes(x=qunif(seq(0,1,length.out=length(quants_sampsplit_split3))), y=quants_sampsplit_split3, col="Sample Splitting"),alpha=1, lwd=1.6)+
  geom_line(aes(x=qunif(seq(0,1,length.out=length(quants_treevalue_split3))), y=quants_treevalue_split3, col="Selective Z-test, Known Var"), alpha=1, lwd=1.1)+
  geom_line(aes(x=qunif(seq(0,1,length.out=length(quants_treevalue_est3))), y=quants_treevalue_est3, col="Selective Z-test, Estimated Var"),alpha=1, lwd=1.1)+
  xlab("") +
  ylab("")+ggtitle("Level 3") +
  geom_abline()+
  theme(plot.title = element_text(size = 10),axis.title=element_text(size = 7.5)) +
  coord_fixed()+theme_bw() + 
  labs(col="method")+
  scale_color_manual(values = c(naivezcol, sampsplitcol, treevalcol, treevalESTcol))

### Note to self!!! The sizes were not working out well automatically. Cropped out white space in preview.
library(patchwork)
p1+p2+p3+plot_layout(guides="collect")
ggsave("~/Dropbox/Tree Values Paper/JMLR Reviewer Responses/type1.png")
