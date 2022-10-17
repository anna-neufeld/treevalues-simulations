setwd("~/Dropbox/Tree Values Paper/Code Reviewer Response")
source("Code to Run Simulations for Paper/permutation_sims.R")
library(intervals)
library(rpart)
library(treevalues)
betas <- rep(seq(0,10, by=1), 200)
n=200
p=10
sigma_y <- 5
maxdepth=3

for (j in 1:length(betas)) {
  print(j)
  filename <- "permutations_3-32-2022.csv"
  
  oneRepFull(n=n,p=p,sigma_y=sigma_y, seed=j, beta=betas[j],
             filename =   filename ,
             minbucket=1, maxdepth=3, XORlev=0.5)
  oneRepFull(n=n,p=p,sigma_y=sigma_y, seed=j, beta=betas[j],
             filename = filename  ,
             minbucket=1, maxdepth=3, XORlev=1)
  oneRepFull(n=n,p=p,sigma_y=sigma_y, seed=j, beta=betas[j],
             filename = filename  ,
             minbucket=1, maxdepth=3, XORlev=2)
}
