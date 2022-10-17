Z.test <- function(dat, alpha, sigmasq) {
  zmult <- qnorm(1-alpha/2)
  confint <- c(mean(dat) - zmult*sqrt(sigmasq/NROW(dat)),
               mean(dat) + zmult*sqrt(sigmasq/NROW(dat)))
  return(list("conf.int" = confint))
}


find_alpha_ideal <- function(n,p,depth, alpha) {
  alphas <-  exp(seq(-10, log(alpha), length.out=50))
  covs <- rep(NA, length(alphas))
  for (k in 1:length(alphas)) {
    alpha_k <- alphas[k]
    print("---------")
    print(alpha_k)
    coverage <- rep(NA, 1000)
    for (trial in 1:1000) {
      set.seed(trial)
      X <- matrix(rnorm(n*p, 0,1), ncol=p)
      y <- rnorm(n, mean=0, sd=1)
      dat <- data.frame(X,y)
      names(dat) <- c(paste("X", 1:p, sep=""), "y")
      form <- paste("y~", paste("X", 1:p, sep="", collapse="+"))
      tree <- rpart::rpart(formula(form), data=dat,model=TRUE, maxdepth=depth, cp=0)
      nodes <- unique(tree$where)
      miniCovs <- rep(NA, length(nodes))
      c <- 1
      for (node in nodes) {
        int <- Z.test(y[tree$where==node], alpha_k, 1)$conf.int[1:2]
        miniCovs[c] <- int[1] < 0 & int[2] > 0
        c <- c+1
      }
      coverage[trial] <- mean(miniCovs)
    }
    covs[k] <- mean(coverage)
    if (covs[k] < 1-alpha-0.1) {break}
    print(covs[k])
  }
  
  ptrue <-  which(covs < 1-alpha)[1]
  ftrue <- (covs[ptrue-1]-1+alpha)/(covs[ptrue-1]-covs[ptrue])
  alpha.prime.true <- (1-ftrue)*alphas[ptrue-1] + ftrue*alphas[ptrue]
  return(alpha.prime.true)
}


alpha.prime.2.1 <- find_alpha_ideal(100,2,1,0.1)
alpha.prime.2.2 <- find_alpha_ideal(100,2,2,0.1)
alpha.prime.2.3 <- find_alpha_ideal(100,2,3,0.1)
alpha.prime.20.1 <- find_alpha_ideal(100,20,1,0.1)
alpha.prime.20.2 <- find_alpha_ideal(100,20,2,0.1)
alpha.prime.20.3 <- find_alpha_ideal(100,20,3,0.1)
alpha.primes <- c(alpha.prime.2.1,
                  alpha.prime.2.2,
                  alpha.prime.2.3,
                  alpha.prime.20.1,
                  alpha.prime.20.2,
                  alpha.prime.20.3)
ps <- c(2,2,2,20,20,20)
depths <- c(1,2,3,1,2,3)


test.cov <- matrix(NA, nrow=2000, ncol=6)
n <- 100
alpha <- 0.1
for (trial in 1:2000) {
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
    
    nodes <- unique(tree$where)
    miniCovs <- rep(NA, length(nodes))
    c <- 1
    for (node in nodes) {
      int <- Z.test(y[tree$where==node], alpha.primes[setting], 1)$conf.int[1:2]
      miniCovs[c] <- int[1] < 0 & int[2] > 0
      c <- c+1
    }
    test.cov[trial, setting] <- mean(miniCovs)
  }
}
  
  
  
print(colMeans(test.cov))
    
    