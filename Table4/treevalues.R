library(treevalues)

test.cov <- matrix(NA, nrow=1000, ncol=6)
n <- 100
alpha <- 0.1
ps <- c(2,2,2,20,20,20)
depths <- c(1,2,3,1,2,3)

for (trial in 1:1000) {
  if (trial%%100==0) {
    save(test.cov, file="treevals_results.RData")
  }
  print(trial)
  set.seed(trial)
  for (setting in 1:6) {
    p <- ps[setting]
    X <- matrix(rnorm(n*p, 0,1), ncol=p)
    y <- rnorm(n, mean=0, sd=1)
    dat <- data.frame(X,y)
    names(dat) <- c(paste("X", 1:p, sep=""), "y")
    form <- paste("y~", paste("X", 1:p, sep="", collapse="+"))
    tree <- rpart::rpart(formula(form), data=dat,model=TRUE, maxdepth=depths[setting], cp=0)
    
    nodes <- row.names(tree$frame)[unique(tree$where)]
    miniCovs <- rep(NA, length(nodes))
    c <- 1
    for (node in nodes) {
      branch <- getBranch(tree, node)
      int <- branchInference(tree,branch,type="reg", computeCI=TRUE, alpha=0.1)$confint
      miniCovs[c] <- int[1] < 0 & int[2] > 0
      c <- c+1
    }
    test.cov[trial, setting] <- mean(miniCovs)
  }
}

setwd("~/Dropbox/Tree Values Paper/Treevalues_Paper_Code/Table4")
load("treevals_results.RData")
mean(is.na(test.cov))
print(colMeans(test.cov, na.rm=TRUE))
