Z.test <- function(dat, alpha, sigmasq) {
  zmult <- qnorm(1-alpha/2)
  confint <- c(mean(dat) - zmult*sqrt(sigmasq/NROW(dat)),
               mean(dat) + zmult*sqrt(sigmasq/NROW(dat)))
  return(list("conf.int" = confint))
}


bootstrap_procedure <- function(tree, B=100, K=10, alpha=0.05, trueMeanVec = NULL, 
                                minExp=-20) {
  gammas <- rep(0,K) # individual coverage probabilities

  y <- tree$model[,1]
  X <- tree$model[,-1]
  n <- NROW(X)
  p <- NCOL(X)

  for (b in 1:B) {
    bootSample <- sample(1:n, size=n, replace=TRUE)
    Xboot <- X[bootSample,]
    yboot <- y[bootSample]

    dat2 <- data.frame(cbind(Xboot, yboot))
    names(dat2) <- c(paste("X", 1:p, sep=""), "y")
    
    tree_b <-  rpart::rpart(tree$call$formula, data=dat2, control=tree$control, model=TRUE)
    root_nodes_b <- unique(tree_b$where)
    
    
    #### I believe that this will make the PREDICT function return WHICH ROW you go in. 
    tree_b$frame$yval <- 1:NROW(tree_b$frame)
    orig_data_root_assignments <- predict(tree_b, newdata=dat)
    Lb <- length(root_nodes_b)
    
    c <- matrix(0, nrow=K, ncol= Lb)
    
    l=0
    alphas <- exp(seq(minExp, log(alpha), length.out=K))
    ### I have no loop over G because I only have an intercept. 
    
    ## Loop over terminal nodes. 
    for (root in root_nodes_b) {
      l = l+1
      B0.boot <- mean(yboot[tree_b$where==root])
      
      ### Shouldn't this be estimated on like everyone who is held OUT of the bootstrap sample?? 
      B0.true <- mean(y[ orig_data_root_assignments==root]) ## not really true lol
     
      #B0.actually.true <- 0
      
      #write(c(B0.boot, B0.true, B0.actually.true, sum(tree_b$where==root)), file="understand_res_inner_boot7", append=TRUE, ncolumns = 4)
      
      
      for (k in 1:K)  {
        interval <- try(Z.test(yboot[tree_b$where==root], alphas[k],1)$conf.int[1:2])
        if (class(interval) != "try-error") {
        if (interval[1] <=  B0.true & interval[2] >= B0.true) {
        c[k,l] <- 1
        }
        }
      }
    }
    
    for (k in 1:K) {
      gammas[k] <- gammas[k] + 1/(Lb)*sum(c[k,])
    }
  }
  
    gammas <- gammas/B
    p <- which(gammas < 1-alpha)[1]
    f <- (gammas[p-1]-1+alpha)/(gammas[p-1]-gammas[p])
    alpha.prime <- (1-f)*alphas[p-1] + f*alphas[p]
    return(alpha.prime)
}



test.cov <- matrix(NA, nrow=1000, ncol=6)
alpha.primes <- matrix(NA, nrow=1000, ncol=6)
n <- 100
alpha <- 0.1
ps <- c(2,2,2,20,20,20)
depths <- c(1,2,3,1,2,3)

for (trial in 1:1000) {
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
    
    alpha.prime <- bootstrap_procedure(tree, K=50, alpha=0.1, minExp=-10)
    alpha.primes[trial, setting] <- alpha.prime
    
    nodes <- unique(tree$where)
    miniCovs <- rep(NA, length(nodes))
    c <- 1
    for (node in nodes) {
      int <- Z.test(y[tree$where==node], alpha.prime, 1)$conf.int[1:2]
      miniCovs[c] <- int[1] < 0 & int[2] > 0
      c <- c+1
    }
    test.cov[trial, setting] <- mean(miniCovs)
  }
}



print(colMeans(test.cov, na.rm=TRUE))
print(colMeans(alpha.primes, na.rm=TRUE))





