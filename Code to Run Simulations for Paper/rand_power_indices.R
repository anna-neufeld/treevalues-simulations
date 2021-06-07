library(rpart)
library(intervals)
library(partykit)
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
  base_tree <- rpart::rpart(
    y~., data=dat, model=TRUE,
    control=rpart.control(maxcompete=0,maxsurrogate=0,
                          xval=0, cp=complex/(sum((y-mean(y))^2)),
                          maxdepth=maxdepth,
                          minbucket=minbucket))



  terminalNodes <- sort(unique(base_tree$where))

  if (length(terminalNodes) > 1) {
    splitList <- getAllBranches(base_tree)
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
      write(paste(c(complex, seed, depth, beta, "condition",p_split,
                    true_signal_split, length(y1), length(y2), sample_signal,
                    rand, i, XORlev, best,'\n'), collapse=" "), file=filename,append=TRUE)

    }
  } else {
    ### If the fitted tree was only a root, it got onone of the true splits!!!
    for (i in 1:4) {
      write(paste(c(complex, seed, NA, beta, "condition",NA,
                    NA, NA,NA, NA, 0, i,XORlev, NA,'\n'), collapse=" "), file=filename,append=TRUE)
    }
  }


  ### Repeat the process, but now on CTree!!!!
  Ctree <- ctree(y~., data=dat, control=ctree_control(minbucket = minbucket,
                                                      maxdepth=maxdepth,
                                                      alpha=alpha))

  preds <- predict(Ctree, type="node")
  internals <- nodeids(Ctree)[-nodeids(Ctree, terminal = TRUE)]
  #terminalnodes <- nodeids(Ctree, terminal = TRUE)
  depths <- table(unlist(sapply(internals, function(x) intersect(internals, nodeids(Ctree, from = x)))))
  names(depths) <- as.character(internals)


  ### Encode all the fitted splits as vecs.
  if (length(internals) > 0) {
    all_wheres <- matrix(nrow=length(internals),ncol=n)
    all_child1 <- rep(0, length(internals))
    all_child2 <- rep(0, length(internals))
    all_nodes <-  rep(0, length(internals))
    j=1

    for (node in internals) {
      children = nodeapply(Ctree, ids =node, function(n) kids_node(n))[[1]]
      all_child1[j] <- children[[1]]$id
      all_child2[j] <- children[[2]]$id
      indices <- c(as.numeric(rownames( Ctree[children[[1]]$id]$dat)),
                   as.numeric(rownames( Ctree[children[[2]]$id]$dat)))
      all_wheres[j,] <- nodeapply(Ctree,id=node,function(n) kidids_node(n,data=dat))[[1]]
      all_wheres[j,][-indices] <- 0
      all_nodes[j] <- node
      j=j+1
    }


    ### For each true split, find best elligible fitted split.
    ### Record inference properties.
    for (i in 1:4) {
      indices <- trueIndices[i,]
      truth <- trueSplits[i,]
      randfulls <- apply(all_wheres,1, function(u) mclust::adjustedRandIndex(u[indices], truth[indices]))
      indices1 <- sapply(all_child1, function(u) as.numeric(rownames( Ctree[u]$dat)),
                         simplify=FALSE)
      indices2 <- sapply(all_child2, function(u) as.numeric(rownames( Ctree[u]$dat)),
                         simplify=FALSE)

      numInvolved <- sapply(1:length(indices1), function(u) length(indices1[[u]])+length(indices2[[u]]))

      elligible <- which(numInvolved < 1.2*trueSize[i] & numInvolved > 0.8*trueSize[i])
      if (length(elligible)==0) {next}
      best <- elligible[which.max(randfulls[elligible])]
      rand <- max(randfulls[elligible])
      node <- all_nodes[[best]]
      depth <- depths[[as.character(node)]]
      child1 <- all_child1[[best]]
      child2 <- all_child2[[best]]
      p_split = as.numeric(nodeapply(Ctree, ids =node, function(n) info_node(n)$p.value))
      y1 <- Ctree[child1]$dat$y
      y2 <- Ctree[child2]$dat$y
      sample_signal <- mean(y1)-mean(y2)
      indice1 <- as.numeric(rownames( Ctree[child1]$dat))
      indice2 <- as.numeric(rownames( Ctree[child2]$dat))
      indices <- c(indice1, indice2)
      true_signal_split <- mean(mu_y[indice1])-mean(mu_y[indice2])

      write(paste(c(complex, seed, depth, beta, "ctree",p_split,
                    true_signal_split, length(y1), length(y2), sample_signal,
                    rand, i, XORlev, best, '\n'), collapse=" "), file=filename,append=TRUE)
    }
  }   else {
    for (i in 1:4) {
      write(paste(c(complex, seed, 1, beta, "ctree",NA,
                    NA, NA,NA, NA, 0, i,XORlev, NA, '\n'), collapse=" "), file=filename,append=TRUE)
    }
  }

  #### SAMPLE SPLITTING
  n1 <- sample(1:n, size=floor(n/2))
  n2 <- setdiff(1:n, n1)
  dat1 <- dat[n1,]
  dat2 <- dat[-n1,]
  split_tree <- rpart::rpart(y~., data=dat1, model=TRUE,
                             control=rpart.control(maxcompete=0,maxsurrogate=0,
                                                   xval=0, cp=complex/(sum((dat1$y-mean(dat1$y))^2)),
                                                   maxdepth=maxdepth,minbucket=minbucket))

  terminalNodes <- sort(unique(split_tree$where))
  split_tree$frame$yval = 1:NROW(split_tree$frame)
  test_predict = predict(split_tree, newdata=dat2)

  if (length(terminalNodes) > 1) {
    splitList <- getAllBranches(split_tree)
    i=1

    ### DEFINITIONS OF DETECTION ARE ON ALL DATA: train + test
    all_wheres <- matrix(ncol=n,nrow=length(splitList)/2)
    all_wheres_test <- matrix(ncol=length(n2), nrow=length(splitList)/2)
    all_splits <- list()
    j=1

    while (i < length(splitList)) {
      splits <- splitList[i][[1]]
      splitText1 <- paste(paste(
        "dat$", splits), collapse=" & ")
      splitText2 <- paste(paste("dat$", splitList[i+1][[1]]), collapse=" & ")
      i <- i+2
      node1 <- eval(parse(text =splitText1))
      node2 <- eval(parse(text =splitText2))*2
      all_wheres[j,] <- node1+node2
      all_wheres_test[j,] <- (node1+node2)[n2]
      all_splits[[j]] <- splits
      j <-j+1
    }

    for (i in 1:4) {
      indices <- trueIndices[i,]
      truth <- trueSplits[i,]
      randfulls <- apply(all_wheres,1, function(u) mclust::adjustedRandIndex(u[indices], truth[indices]))
      elligible <- which(rowSums(all_wheres != 0) < 1.2*trueSize[i] & rowSums(all_wheres != 0) > 0.8*trueSize[i])
      if (length(elligible)==0) {next}
      best <- elligible[which.max(randfulls[elligible])]
      rand <- max(randfulls[elligible])
      splits <- all_splits[[best]]
      where <- all_wheres[best,]
      where_test <- all_wheres_test[best,]

      ##### DO INFERENCE ON TEST SET
      y1 <- dat2[where_test==1,]$y
      y2 <- dat2[where_test==2,]$y

      nu_full <- (where==1)/sum(where==1) - (where==2)/sum(where==2)
      true_signal_split <- t(nu_full)%*%mu_y
      sample_signal_split <- mean(y1)-mean(y2)
      se <- sqrt(sigma_y^2/length(y1) + sigma_y^2/length(y2))
      p_split <- (1-pnorm(abs((mean(y1)-mean(y2))/se)))*2
      depth = length(splits)


      write(paste(c(complex, seed, depth, beta, "splitting",p_split,
                    true_signal_split, length(y1), length(y2), sample_signal_split,
                    rand, i, XORlev, best, '\n'), collapse=" "), file=filename,append=TRUE)
    }
  } else {
    for (i in 1:4) {
      write(paste(c(complex, seed, NA, beta, "splitting", NA,
                    NA, NA,NA, NA, 0, i,XORlev, NA, '\n'), collapse=" "), file=filename,append=TRUE)
    }
  }
}


















