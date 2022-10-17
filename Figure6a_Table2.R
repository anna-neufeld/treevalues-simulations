library(tidyverse)
library(gridExtra)

naivezcol <- "#7CAE00"
sampsplitcol <- "#00BFC4"
treevalcol <- "#F8766D"
ctreecol <- "black"

## Simulations for this are in file "non_null_CI_sims.R" and are run from "non_null_sims_RUN.R". These are the only set of sims that are really
## slow because CIs are way slower than pvals. Totally should have been run on cluster but, whoops, just did on computer.
setwd("~/Dropbox/Tree Values Paper/Code : Other for Creating All Figures")
res <- read.csv("non_null_CI_main-1-18-2021", header=FALSE, sep=" ")
names(res) <- c("beta", "XORLev", "cp", "seed", "method", "type", "depth", "pval", "lower", "upper", "n1", "n2", "sampsig","truesig", "truesigtest", "sizeS", "sizeNum")
res <- res %>% mutate(correctCI = lower < truesig & upper > truesig, correctCItest = lower < truesigtest & upper > truesigtest,
                      length = upper - lower, invalid = n1*n2==0)
res <- res %>% filter(XORLev != 0, XORLev !=3)


###  TABLE 2 !!!
res %>% group_by(method, type, depth) %>% summarize("Full Coverage" = round(mean(correctCI, na.rm=TRUE),3),"Test Coverage" = round(mean(correctCItest, na.rm=TRUE),3))



### FIGURE 6
groupRes <- res %>%
  filter(depth < 4,XORLev %in% c(0.5,1,2)) %>%
  group_by(method, beta, XORLev, depth,type) %>% summarize(length=median(length,na.rm=TRUE))
groupRes <- groupRes %>% mutate(XORlev2 = paste0("a=", XORLev),
                                depth2 = paste0("Level ", depth),
                                method=ifelse(method=="splitting", "Sample splitting", method))

groupRes$method2 <- "(??), reg"
groupRes$method2[groupRes$method=="Tree-Values" & groupRes$type=="split"] <-
  "(10), sib"
groupRes$method2[groupRes$method=="Sample splitting" & groupRes$type=="split"] <-
  "Sample splitting, sib"
groupRes$method2[groupRes$method=="Sample splitting" & groupRes$type=="child"] <-
  "Sample splitting, reg"

p4 <- ggplot(data=groupRes %>% filter(method=="Tree-Values"), aes(x=beta, y=length, col=method2, lty=method2))+ geom_smooth()+facet_grid(cols=vars(XORlev2), rows=vars(depth2))+scale_y_log10()+ylab("Median Confidence Interval Width")+xlab("b")+
  scale_colour_manual(
    name="Method",
    values = c(treevalcol,treevalcol,sampsplitcol, sampsplitcol),
    labels=c(
      expression(paste("Selective Z-interval in (23),", nu[reg]^T~mu, "        ")),
      expression(paste("Selective Z-interval in (12),", nu[sib]^T~mu, "        "))))+
  theme_bw()+ scale_linetype_manual(
    name = "Method", values = c(1,3,1,3),
    labels=c(
      expression(paste("Selective Z-interval in (29),", nu[reg]^T~mu, "        ")),
      expression(paste("Selective Z-interval in (11),", nu[sib]^T~mu, "        "))))+
  theme(axis.title=element_text(size=12), title=element_text(size=12), axis.text=element_text(size=12))


p4log <- ggplot(data=groupRes %>% filter(method=="Tree-Values",type=="child"), 
                aes(x=beta, y=length, col=method2, lty=depth2))+ 
  geom_smooth(se=FALSE, lwd=1)+facet_grid(cols=vars(XORlev2))+ylab("Median Confidence Interval Width")+xlab("b")+
  theme_bw()+
  scale_linetype_manual(values=c(1,2,3))+theme(axis.title=element_text(size=12), title=element_text(size=12), axis.text=element_text(size=12))

p4log + guides(col=FALSE, lty=FALSE)
ggsave(filename="~/Dropbox/Tree Values Paper/JMLR_resubmit_Oct_2022/Figures/CI_width_nolog_nolev.png")
