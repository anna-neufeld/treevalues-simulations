library(tidyverse)
library(gridExtra)


#naivezcol <- "#7CAE00"
#sampsplitcol <- "#00BFC4"
treevalcol <- "#F8766D"
sigmaSSEcol <- "navy"
sigmaConscol <- "darkgray"
#treevalESTcol <- "red"
#ctreecol <- "black"

setwd("~/Dropbox/Tree Values Paper/Code Reviewer Response/unknown variance/")
nullRes <- read.csv("null_res_var_est.csv", sep=" ", header=FALSE)
names(nullRes) <- c("cp", "seed", "method", "type", "depth", "pval", "n1", "n2", "sampsig", "NA")

nullRes <- nullRes %>% filter(method != "mad")
nullRes$method = ordered(nullRes$method, levels=c("cons", 
                                                  "SSE", "known", "med_yc"))

nullRes$level = paste("Level", nullRes$depth)

ggplot(data=nullRes %>% filter(method != 'med_yc'), aes(sample=pval, col=method))+
  geom_qq(distribution="qunif", alpha=0.8)+
  facet_grid(cols=vars(level))+
  geom_abline(slope=1, intercept=0)+coord_fixed()+theme_bw()+
  xlab("Unif(0,1) Quantiles")+ylab("Empirical Quantiles")+
  labs(col="")+
  scale_color_manual(labels=c(expression(hat(sigma)[cons]), expression(hat(sigma)[SSE]),
                                expression(paste("True ", sigma))),
                     values=c(sigmaConscol, sigmaSSEcol, treevalcol))+
  theme(legend.text.align = 0)
  #theme(strip.text = element_text(size = 20))
ggsave("~/Dropbox/Tree Values Paper/JMLR Reviewer Responses/global_null_fig.png",
       width=8, height=4)
ggsave("~/Dropbox/Tree Values Paper/JMLR_submit_v3/Figures/global_null_fig_est.png",
       width=8, height=4)
