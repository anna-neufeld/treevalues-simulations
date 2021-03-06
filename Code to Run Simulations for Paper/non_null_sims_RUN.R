setwd("~/treevalues")
devtools::load_all()


setwd("~/REPRODUCE FIGURES TO MAKE SURE THEY STILL WORK/")
source("Code to Run Simulations for Paper/non_null_CI_sims.R")



library(intervals)
library(rpart)
betas <- rep(seq(0,10, by=1), 500)
n=200
p=10
sigma_y <- 5
maxdepth=3

for (j in 1:length(betas)) {
  print(j)
  filename <- "non_null_CI_main-5-27-2021"

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
