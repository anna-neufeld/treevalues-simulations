library(rpart)
library(intervals)
#library(partykit)
library(rpart.utils)


### Compares CTree, Sample Splitting, and Selective Z Tests in terms of power.
### Also needs to store adjusted rand index to figure out which splits in the tree match the TRUE splits.
full_Comparison_NEW <-  function(n,p,sigma_y, seed, complex = 0.01*n*sigma_y^2, alpha = 0.05, beta=0, filename="test.txt",
                                 minbucket=7,maxdepth=5, XORlev=1)  {


  set.seed(seed)


  ### Generate data accordng to provided parameters.
  ### beta is called b in the paper.
  ### XORlev is called a in the paper.
  X <- MASS::mvrnorm(n, rep(0,p), diag(rep(1,p)))
  mu_y <- beta*I(X[,1] < 0) +
    beta*XORlev*(I(X[,1] < 0 & X[,2] > 0)) +
    beta*I(X[,3] < 0 & X[,2] < 0 & X[,1] < 0)+
    beta*I(X[,3] > 0 & X[,2] > 0 & X[,1] < 0)
  y <- rnorm(n,mu_y,sigma_y)
  dat <- data.frame(y=y,X=X)
  nameX <- sapply(1:p, function(u) paste0("X",u))
  names(dat) = c("y", nameX)

  ### Get vectors that encode the 4 TRUE splits in the data.
  split1true <- I(X[,1] < 0)
  split2true <- split1true + I(X[,1] < 0  & X[,2]>0)
  split3true <- split2true + 2*I(X[,3] > 0 & X[,2] < 0 & X[,1] < 0)
  split3true[split3true==2] <- 0
  split4true <- split2true + I(X[,3] > 0 & X[,2] > 0 & X[,1] < 0) + 2*I(X[,3] < 0 & X[,2] > 0 & X[,1] < 0)
  split4true[split4true < 3] <- 0


  trueSplits <- rbind(
    split1true,
    split2true,
    split3true,
    split4true
  )

  ### Encodes which observations are actually INVOLVE in the true splits.
  trueIndices <- rbind(
    rep(TRUE,n),
    split2true > 0,
    split3true != 0,
    split4true != 0
  )

  trueSize <- rowSums(trueIndices != 0)

  #### TREE VALUES CALCS ########################################
  base_tree <- rpart::rpart(y~., data=dat, model=TRUE,
    control=rpart.control(maxcompete=0,maxsurrogate=0,
                          xval=0, cp=complex/(sum((y-mean(y))^2)),
                          maxdepth=maxdepth,
                          minbucket=minbucket))



  terminalNodes <- sort(unique(base_tree$where))

  if (length(terminalNodes) > 1) {
    splitList <- getBranch(base_tree)
    i=1

    all_wheres <- matrix(ncol=n,nrow=length(splitList)/2)
    all_splits <- list()
    j=1

    ### Store each split in the fitted tree as a vector!!!!
    while (i < length(splitList)) {
      splits <- splitList[i][[1]]

      splitText1 <- paste(paste("dat$", splits), collapse=" & ")
      splitText2 <- paste(paste("dat$", splitList[i+1][[1]]), collapse=" & ")
      i <- i+2
      node1 <- eval(parse(text =splitText1))
      node2 <- eval(parse(text =splitText2))*2
      all_wheres[j,] <- node1+node2
      all_splits[[j]] <- splits
      j <-j+1
    }

    ### For each TRUE split, figure out which FITTED split has the highest RAND INDEX match.
    ### A split is only elligible to be a match if approximately the right number of observations are involved.
    ### This is to avoid settings where the rand index looks artificially good because the fitted split like only involves 2 observations.
    for (i in 1:4) {
      indices <- trueIndices[i,]
      truth <- trueSplits[i,]
      randfulls <- apply(all_wheres,1, function(u) mclust::adjustedRandIndex(u[indices], truth[indices]))
      elligible <- which(rowSums(all_wheres != 0) < 1.2*trueSize[i] & rowSums(all_wheres != 0) > 0.8*trueSize[i])
      if (length(elligible)==0) {next}

      ### For the best elligible split, record the pvalue for the split!! And also record the rand!
      ### and a few other things like true sig, sample size, etc. just for fun,
      best <- elligible[which.max(randfulls[elligible])]
      rand <- max(randfulls[elligible])
      splits <- all_splits[[best]]
      where <- all_wheres[best,]
      nu <- (where==1)/sum(where==1) - (where==2)/sum(where==2)
      true_signal_split <- t(nu)%*%mu_y
      y1 <- y[where==1]
      y2 <- y[where==2]
      sample_signal <- t(nu)%*%y
      depth = length(splits)
      phi_bounds_split <- getInterval(base_tree, nu,splits)
      p_split <- correctPVal(phi_bounds_split, nu, y, sigma_y)
      
      p_split_cons <- correctPVal(phi_bounds_split, nu, y, sd(y))
      p_split_SSE <- correctPVal(phi_bounds_split, nu, y, sqrt(sum((y-predict(base_tree))^2)/(n-length(unique(base_tree$where)))))
      p_split_med_yc <- correctPVal(phi_bounds_split, nu, y, sqrt(median((y - median(y))^2)/qchisq(0.5, df=1)))
      p_split_mad <- correctPVal(phi_bounds_split, nu, y, median(abs(y - median(y))))
      
      
      write(paste(c(complex, seed, depth, beta, "condition",p_split, p_split_cons, p_split_SSE, p_split_med_yc, 
                    p_split_mad,
                    true_signal_split, length(y1), length(y2), sample_signal,
                    rand, i, XORlev, best,'\n'), collapse=" "), file=filename,append=TRUE)
    

    }
  } else {
    ### If the fitted tree was only a root, it got onone of the true splits!!!
    for (i in 1:4) {
      write(paste(c(complex, seed, NA, beta, "condition",NA,NA,NA,NA,NA
                    NA, NA,NA, NA, 0, i,XORlev, NA,'\n'), collapse=" "), file=filename,append=TRUE)
    }
  }
}


















