setwd("~/Dropbox/Tree Values Paper/Code : Other for Creating All Figures/")
library(rpart)
library(partykit)
library(treevalues)

## While the blsdata dataset is now built into the treevalues package, it was originally from the
## visTree github page.


## This code comes firectly from the VisTree Github Page
covariates <-blsdata %>% dplyr::select(hunger, disinhibition, resteating, rrvfood, liking, wanting, kcal24h0)
covariates<-data.frame(scale(covariates))
p2<-partykit::ctree(kcal24h0~hunger+ disinhibition+resteating+rrvfood+ liking+wanting, data = covariates, control = ctree_control(mincriterion = 0.95))
p4<-rpart(kcal24h0~., model = TRUE, data = covariates,control = rpart.control(cp = 0.0185, maxcompete = 0, maxsurrogate = 0))


p4$frame$yval[1] <- 0 ## Avoid printing 10^{-18} when I for sure know that it is supposed to be a 0

# This makes a throwaway rpart object with the same size as the fitted CTree tree p2.
# It is a hacky way to manually turn our CTree into an rpart object.
temp <-rpart(kcal24h0~., model = TRUE, data = covariates,control = rpart.control(maxdepth=3,cp=0.03001, maxcompete = 0, maxsurrogate = 0))
temp <- inferenceFrame(temp)
temp$frame$pval <- c(NA, "< 0.001", "< 0.001", "< 0.001", "= 0.016", "= 0.016", "<0.01")

## Now manually (using p2 from above as a guide)-- turn this rpart into the CTree!!
## Make the splits say the correct thing for labeling splits!!
temp$frame$var[2] <- "liking"
row.names(temp$splits)[2] <- "liking"
temp$splits[1,4] <- 1.849
temp$splits[2,4] <- -0.277
temp$splits[3,4] <- -1.261

#### Figure out all of the fitted means from the CTree object so that they can be pasted in.
#### The point of this is just to be able to plot the CTree using our own plotting function.
mean(covariates$kcal24h0[covariates$hunger <= 1.693])
mean(covariates$kcal24h0[covariates$hunger > 1.693])
mean(covariates$kcal24h0[covariates$hunger <= 1.693 & covariates$liking <= -0.277])
mean(covariates$kcal24h0[covariates$hunger <= 1.693 & covariates$liking > -0.277])
mean(covariates$kcal24h0[covariates$hunger <= 1.693 & covariates$liking > -0.277 & covariates$rrvfood <= -1.261])
mean(covariates$kcal24h0[covariates$hunger <= 1.693 & covariates$liking > -0.277 & covariates$rrvfood > -1.261])

temp$frame$yval<- c(0, -0.09118071,-0.43876,0.1414685,-0.5060532,0.234805,1.38073645)


#### FINALLY READY TO MAKE THE THREE PANELS OF THE PLOT.

##### NOTE TO SELF::: TO MAKE THE TREE-PLOTS SQUATTER IN THE PAPER, I MANUALLY CROPPED THE BRANCH LINES SHORTER IN PREVIEW
n <- NROW(blsdata)
nT <- length(unique(p4$where))
SSE <- sqrt(1/(n-nT)*sum((covariates$kcal24h0 - predict(p4))^2))
#sd(blsdata$kcal24h0)

png(filename = "~/Dropbox/Tree Values Paper/JMLR Reviewer Responses/bls_tree_SSE.png", width=3, height=5,
    units="in", res=400)
treeval.plot(p4, sigma_y = SSE, inferenceType=4,alpha=0.05, digits=2, printn = FALSE, nn=FALSE)
dev.off()

png(filename = "~/Dropbox/Tree Values Paper/JMLR Reviewer Responses/bls_tree_cons.png", width=3, height=5,
    units="in", res=400)
treeval.plot(p4,inferenceType=4,alpha=0.05, digits=2, printn = FALSE, nn=FALSE)
dev.off()




