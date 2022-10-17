library(rpart)
library(intervals)
library(partykit)
library(rpart.utils)
library(treevalues)

#### Im gonna shut off pruning and use a max depth of 3
#### The pruning calculations are like truly dumb
#### Although I *do* want someone to believe I could do them if I wanted to, so was this a bad move??


oneRepFull <-  function(n,p,sigma_y, seed, complex = 0, alpha = 0.05, beta=0, filename="test.txt",
                        minbucket=1,maxdepth=3, XORlev=1) {

  set.seed(seed)
  X <- MASS::mvrnorm(n, rep(0,p), diag(rep(1,p)))

  mu_y <- beta*I(X[,1] < 0) +
    beta*XORlev*(I(X[,1] < 0 & X[,2] > 0)) +
    beta*I(X[,3] < 0 & X[,2] < 0 & X[,1] < 0)+
    beta*I(X[,3] > 0 & X[,2] > 0 & X[,1] < 0)
  
  y <- rnorm(n,mu_y,sigma_y)
  dat <- data.frame(y=y,X=X)
  nameX <- sapply(1:p, function(u) paste0("X",u))
  names(dat) = c("y", nameX)

  ### Build rpart tree.
  base_tree <- rpart::rpart(y~., data=dat,
                          control=rpart.control(minbucket=minbucket,cp=0, maxcompete=0,maxsurrogate=0, maxdepth=3), model=TRUE)


  if (length(unique(base_tree$where)) > 1) {
    #### Go through EVERY SPLIT and EVERY NODE
    splitList <- getBranch(base_tree)
    lengths <- lapply(splitList, length)
    indices <- which(lengths==3)
    if (length(indices)>0) {
      splitList <- splitList[indices]
      for (branch in splitList) {
        ########## COMPUTE A PVAL FOR EACH PERMUTATION
        branch2 <- paste("dat$", branch)
        splitText <- paste(branch2, collapse=" & ")
        node1 <- eval(parse(text =splitText))
        nu <- (node1==1)/sum(node1==1)
        sample_signal <- t(nu)%*%y
        true_signal <- t(nu)%*%mu_y
        list1 <- combinat::permn(branch)
      
        all_bounds <- list()
        counter <- 1
        for (perm in list1) {
          all_bounds[[counter]] <- getInterval(base_tree, nu, perm,sib=FALSE, grow=TRUE)
          counter <- counter+1
        }
        
        lengths <- as.numeric(lapply(all_bounds, length))
        all_bounds <- all_bounds[lengths!=0]
        
        powerSet <- rje::powerSet(all_bounds)
        powerSet <- powerSet[2:length(powerSet)]
        
        merged_bounds <- all_bounds
        for (i in 2:length(all_bounds)) {
          merged_bounds[[i]] <- interval_union(merged_bounds[[i-1]], merged_bounds[[i]])
          closed(merged_bounds[[i]]) <- c(TRUE,TRUE)
          merged_bounds[[i]] <- reduce(merged_bounds[[i]])
        }
        
        lengths0 <- as.numeric(lapply(all_bounds, function(u) sum(size(u))))
        lengths <- as.numeric(lapply(merged_bounds, function(u) sum(size(u))))
        pvals <- rep(NA, length(merged_bounds))
        CIs <- matrix(NA, nrow=length(merged_bounds), ncol=2)
        for (i in 1:length(merged_bounds)) {
          int <-  merged_bounds[[i]] 
          if (lengths0[i] != 0) {
            pvals[i] <- correctPVal(int,nu,y,sigma_y)
            CIs[i,] <- computeCI(nu,y,sigma_y,int)
          } else {
            pvals[i] <-pvals[i-1]
            CIs[i,] <- CIs[i-1,]
          }
        }
        
        res <- c(beta, XORlev, seed, 3, sample_signal, true_signal,
                       lengths0[1], lengths[1], pvals[1], CIs[1,],
                       lengths0[2], lengths[2], pvals[2],CIs[2,],
                       lengths0[3], lengths[3], pvals[3],CIs[3,],
                       lengths0[4], lengths[4], pvals[4],CIs[4,],
                       lengths0[5], lengths[5], pvals[5],CIs[5,],
                       lengths0[6], lengths[6], pvals[6],CIs[6,])
        write(paste(res, collapse=" "), file=filename, append=TRUE)
      }
    }
  }
}
            
        
        
   
        
        
    