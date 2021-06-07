setwd("~/treevalues")
devtools::load_all()
library(intervals)
library(rpart)

setwd("~/REPRODUCE FIGURES TO MAKE SURE THEY STILL WORK/")
source("Code to Run Simulations for Paper/rand_power_indices.R")

betas <- rep(seq(0,10, by=1), 500)

n=200
p=10
sigma_y <- 5
minbucket=1
maxdepth=3
alpha=200

for (j in 1:length(betas)) {
  print(j)
  filename <- "Full_Rand_Comps_5-31-2021-maxdepth3.csv"

  full_Comparison_NEW(n=n,p=p,sigma_y=sigma_y, seed=j, complex=alpha, beta=betas[j],
                      filename =   filename ,
                      minbucket=1, maxdepth=3, XORlev=0.5)
  full_Comparison_NEW(n=n,p=p,sigma_y=sigma_y, seed=j, complex=alpha, beta=betas[j],
                     filename = filename  ,
                      minbucket=1, maxdepth=3, XORlev=1)
  full_Comparison_NEW(n=n,p=p,sigma_y=sigma_y, seed=j, complex=alpha, beta=betas[j],
                      filename = filename  ,
                      minbucket=1, maxdepth=3, XORlev=2)
}
