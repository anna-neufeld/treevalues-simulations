# Simulations are run in the file "Code to Run Simulations for Paper / Null_testing_sims.R"/


library(tidyverse)
library(gridExtra)

naivezcol <- "#7CAE00"
sampsplitcol <- "#00BFC4"
treevalcol <- "#F8766D"
ctreecol <- "black"

setwd("~/Dropbox/Tree Values Paper/Code : Other for Creating All Figures")
nullRes <- read.csv("null_res_1-18-2021-5000.csv", sep=" ", header=FALSE)
names(nullRes) <- c("cp", "seed", "method", "type", "depth", "pval", "n1", "n2", "sampsig", "NA")


quants_treevalue_split1 <- sort((nullRes %>% filter(depth==1,method=="Tree-Values", type=="split"))$pval)
quants_pvalue_split1 <- sort((nullRes %>% filter(depth==1,method=="NaiveZ",type=="split"))$pval)
quants_sampsplit_split1 <- sort((nullRes %>% filter(depth==1,method=="splitting", type=="split"))$pval)

quants_treevalue_split2 <- sort((nullRes %>% filter(depth==2,method=="Tree-Values", type=="split"))$pval)
quants_pvalue_split2 <- sort((nullRes %>% filter(depth==2,method=="NaiveZ",type=="split"))$pval)
quants_sampsplit_split2 <- sort((nullRes %>% filter(depth==2,method=="splitting", type=="split"))$pval)

quants_treevalue_split3 <- sort((nullRes %>% filter(depth==3,method=="Tree-Values", type=="split"))$pval)
quants_pvalue_split3 <- sort((nullRes %>% filter(depth==3,method=="NaiveZ",type=="split"))$pval)
quants_sampsplit_split3 <- sort((nullRes %>% filter(depth==3,method=="splitting", type=="split"))$pval)

p1 <- ggplot(data=NULL) +
  geom_line(aes(x=qunif(seq(0,1,length.out=length(quants_pvalue_split1))), y=quants_pvalue_split1, col="Naive Z-test"), alpha=1, lwd=1.4) +
  geom_line(aes(x=qunif(seq(0,1,length.out=length(quants_sampsplit_split1))), y=quants_sampsplit_split1, col="Sample Splitting"),alpha=1, lwd=1.6)+
  geom_line(aes(x=qunif(seq(0,1,length.out=length(quants_treevalue_split1))), y=quants_treevalue_split1, col="Selective Z-test in (6)"),alpha=1, lwd=1.1)+
  xlab("") +
  ylab("Empirical Quantiles")+ggtitle("Level 1") +
  theme(plot.title = element_text(size = 10),axis.title=element_text(size = 7.5)) +
  coord_fixed()+theme_bw() + labs(col="Method")+guides(col=FALSE)+
  scale_color_manual(values = c(naivezcol, sampsplitcol, treevalcol))

p2 <- ggplot(data=NULL) +
  geom_line(aes(x=qunif(seq(0,1,length.out=length(quants_pvalue_split2))), y=quants_pvalue_split2, col="Naive Z-test"), alpha=1, lwd=1.4) +
  geom_line(aes(x=qunif(seq(0,1,length.out=length(quants_sampsplit_split2))), y=quants_sampsplit_split2, col="Sample Splitting"),alpha=1, lwd=1.6)+
  geom_line(aes(x=qunif(seq(0,1,length.out=length(quants_treevalue_split2))), y=quants_treevalue_split2, col="Selective Z-test in (6)"), alpha=1, lwd=1.1) +
  xlab("U(0,1) Quantiles") +
  ylab("")+ggtitle("Level 2") +
  theme(plot.title = element_text(size = 10),axis.title=element_text(size = 7.5)) +
  coord_fixed()+theme_bw() + labs(col="Method")+guides(col=FALSE)+
  scale_color_manual(values = c(naivezcol, sampsplitcol, treevalcol))

p3 <- ggplot(data=NULL)+
  geom_line(aes(x=qunif(seq(0,1,length.out=length(quants_pvalue_split3))), y=quants_pvalue_split3, col="Naive Z-test"), alpha=1, lwd=1.4) +
  geom_line(aes(x=qunif(seq(0,1,length.out=length(quants_sampsplit_split3))), y=quants_sampsplit_split3, col="Sample Splitting"),alpha=1, lwd=1.6)+
  geom_line(aes(x=qunif(seq(0,1,length.out=length(quants_treevalue_split3))), y=quants_treevalue_split3, col="Selective Z-test in (6)"), alpha=1, lwd=1.1)+
  xlab("") +
  ylab("")+ggtitle("Level 3") +
  theme(plot.title = element_text(size = 10),axis.title=element_text(size = 7.5)) +
  coord_fixed()+theme_bw() + guides(col=FALSE)+
  scale_color_manual(values = c(naivezcol, sampsplitcol, treevalcol))

### Note to self!!! The sizes were not working out well automatically. Cropped out white space in preview.
grid.arrange(p1,p2, p3, widths=c(3,3,3))
#ggsave(filename="~/Dropbox/Tree Values Paper/newformat_v35/Figures/Type1Error_nolegend.png",
#       grid.arrange(p1,p2, p3, widths=c(3,3,3)))
