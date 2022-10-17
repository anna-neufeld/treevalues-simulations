library(treevalues)


setwd("~/Dropbox/Tree Values Paper/Code Reviewer Response/unknown variance")
source("ADD_RAND_TO_CI.R")



library(intervals)
library(rpart)
betas <- rep(seq(0,10, by=1), 500)
n=200
p=10
sigma_y <- 5
maxdepth=3

for (j in 1:length(betas)) {
  print(j)
  filename <- "non_null_var_est_10-10-2022.csv"

  oneRepFull(n=n,p=p,sigma_y=sigma_y, seed=j, complex=200, beta=betas[j],
                      filename =   filename ,
                      minbucket=1, maxdepth=3, XORlev=0.5)
  oneRepFull(n=n,p=p,sigma_y=sigma_y, seed=j, complex=200, beta=betas[j],
                      filename = filename  ,
                      minbucket=1, maxdepth=3, XORlev=1)
  oneRepFull(n=n,p=p,sigma_y=sigma_y, seed=j, complex=200, beta=betas[j],
             filename = filename  ,
             minbucket=1, maxdepth=3, XORlev=2)
}
