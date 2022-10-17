library(rpart)
library(intervals)
library(partykit)
library(rpart.utils)
library(treevalues)


### Runs a giant repetition for the confidence interval coverage and width experiments.
oneRepFull <-  function(n,p,sigma_y, seed, complex = 0.01*n*sigma_y^2, alpha = 0.05, beta=0, filename="test.txt",
                        minbucket=1,maxdepth=5, XORlev=1)
{

  ### Generate dataset according to true mean model using parameters beta, sigma_y, n, and p.
  ### Note that "beta" is called "b" in the paper.
  ### "XORlev" is called "a" in the paper.
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



  ### This is a conversion between "complex", which is the user-provided complexity parameter (called lambda in paper),
  ### and rpart's "cp" parameter, which is a scaled version of lambda.
  cp = complex/sum((y-mean(y))^2)

  ### Build rpart tree.
  base_tree_1 <- rpart::rpart(y~., data=dat,
                            control=rpart.control(minbucket=minbucket,cp=cp, maxcompete=0,maxsurrogate=0, maxdepth=1), model=TRUE)
  base_tree_2 <- rpart::rpart(y~., data=dat,
                              control=rpart.control(minbucket=minbucket,cp=cp, maxcompete=0,maxsurrogate=0, maxdepth=2), model=TRUE)
  
  base_tree <- rpart::rpart(y~., data=dat,
    control=rpart.control(minbucket=minbucket,cp=cp, maxcompete=0,maxsurrogate=0, maxdepth=maxdepth), model=TRUE)

  full_rand_1 <- mclust::adjustedRandIndex(beta*I(X[,1] < 0) , base_tree_1$where)
  full_rand_2 <- mclust::adjustedRandIndex( beta*I(X[,1] < 0) + beta*XORlev*(I(X[,1] < 0 & X[,2] > 0)),
                                            base_tree_2$where)
  full_rand_3 <- mclust::adjustedRandIndex(mu_y, base_tree$where)

  ### If it build a non-trivial rpart tree,
  if (length(unique(base_tree$where)) > 1) {
  #### Go through EVERY SPLIT and EVERY NODE
  splitList <- getBranch(base_tree)
  j=1
  while (j < length(splitList)) {
      splits <- splitList[j][[1]]
      splitText1 <- paste(paste("dat$", splits), collapse=" & ")
      splitText2 <- paste(paste("dat$", splitList[j+1][[1]]), collapse=" & ")
      node1 <- eval(parse(text =splitText1))
      node2 <- eval(parse(text =splitText2))*2
      where <- node1+node2

      y1 <- y[where==1]
      y2 <- y[where==2]

      ### Building nu vectors for both the split and for each individual child.
      nu <- (where==1)/sum(where==1) - (where==2)/sum(where==2)
      nu1 <- (where==1)/sum(where==1)
      nu2<- (where==2)/sum(where==2)
      sample_signal_split <- as.numeric(t(nu)%*%y)
      sample_signal_1 <- as.numeric(t(nu1)%*%y)
      sample_signal_2 <- as.numeric(t(nu2)%*%y)

      true_signal_split <- as.numeric(t(nu)%*%mu_y)
      true_signal_1 <- as.numeric(t(nu1)%*%mu_y)
      true_signal_2 <- as.numeric(t(nu2)%*%mu_y)

      ### SELECTIVE Z INFERENCE. On the split and on the two individual children.
      phi_bounds_split <- getInterval(base_tree, nu,splits)
      phi_bounds_1 <- getInterval(base_tree, nu1,splits)
      phi_bounds_2 <- getInterval(base_tree, nu2,splitList[[j+1]])
      depth = length(splits)
      p_split <- correctPVal(phi_bounds_split, nu, y, sigma_y)
      p_1 <- correctPVal(phi_bounds_1, nu1, y, sigma_y)
      p_2 <- correctPVal(phi_bounds_2, nu2, y, sigma_y)
      
      #p_split_cons <- correctPVal(phi_bounds_split, nu, y, sd(y))
      #p_split_SSE <- correctPVal(phi_bounds_split, nu, y, sqrt(sum((y-predict(base_tree))^2/(n-length(unique(base_tree$where))))))
      #p_split_med_yc <- correctPVal(phi_bounds_split, nu, y, sqrt(median((y - median(y))^2)/qchisq(0.5, df=1)))
      #p_split_mad <- correctPVal(phi_bounds_split, nu, y, median(abs(y - median(y))))
      
      #p_1_cons <- correctPVal(phi_bounds_1, nu1, y, sd(y))
      #p_1_SSE <- correctPVal(phi_bounds_1, nu1, y, sqrt(sum((y-predict(base_tree))^2/(n-length(unique(base_tree$where))))))
      #p_1_med_yc <- correctPVal(phi_bounds_1, nu1, y, sqrt(median((y - median(y))^2)/qchisq(0.5, df=1)))
      #p_1_mad <- correctPVal(phi_bounds_1, nu1, y, median(abs(y - median(y))))
      
      #p_2_cons <- correctPVal(phi_bounds_2, nu2, y, sd(y))
      #p_2_SSE <- correctPVal(phi_bounds_2, nu2, y, sqrt(sum((y-predict(base_tree))^2/(n-length(unique(base_tree$where))))))
      #p_2_med_yc <- correctPVal(phi_bounds_2, nu2, y, sqrt(median((y - median(y))^2)/qchisq(0.5, df=1)))
      #p_2_mad <- correctPVal(phi_bounds_2, nu2, y, median(abs(y - median(y))))
                                 
                                 
      
      CI_split <- computeCI(nu, y, sigma_y,phi_bounds_split, 0.05)
      CI_1 <- computeCI(nu1, y, sigma_y,phi_bounds_1, 0.05)
      CI_2 <- computeCI(nu2, y, sigma_y,phi_bounds_2, 0.05)
      
      #CI_split_cons <- computeCI(nu, y, sd(y),phi_bounds_split, 0.05)
      #CI_1_cons <- computeCI(nu1, y, sd(y),phi_bounds_1, 0.05)
      #CI_2_cons <- computeCI(nu2, y, sd(y),phi_bounds_2, 0.05)
      
      #CI_split_SSE <- computeCI(nu, y, sqrt(sum((y-predict(base_tree))^2)/(n-length(unique(base_tree$where)))),phi_bounds_split, 0.05)
      #CI_1_SSE <- computeCI(nu1, y, sqrt(sum((y-predict(base_tree))^2)/(n-length(unique(base_tree$where)))),phi_bounds_1, 0.05)
      #CI_2_SSE <- computeCI(nu2, y, sqrt(sum((y-predict(base_tree))^2)/(n-length(unique(base_tree$where)))),phi_bounds_2, 0.05)
      
      #CI_split_mad <- computeCI(nu, y, median(abs(y - median(y))),phi_bounds_split, 0.05)
      #CI_1_mad <- computeCI(nu1, y, median(abs(y - median(y))),phi_bounds_1, 0.05)
      #CI_2_mad <- computeCI(nu2, y, median(abs(y - median(y))),phi_bounds_2, 0.05)
      
      #CI_split_med_yc  <- computeCI(nu, y, sqrt(median((y - median(y))^2)/qchisq(0.5, df=1)),phi_bounds_split, 0.05)
      #CI_1_med_yc  <- computeCI(nu1, y, sqrt(median((y - median(y))^2)/qchisq(0.5, df=1)),phi_bounds_1, 0.05)
      #CI_2_med_yc  <- computeCI(nu2, y, sqrt(median((y - median(y))^2)/qchisq(0.5, df=1)),phi_bounds_2, 0.05)

      ### Saving stuff about size of set S. Used to be used for plots.
      ### Didn't end up in any plots in final paper, but still saved bc interesting to look at.
      numerator_split <- suppressWarnings(interval_intersection(interval_complement(phi_bounds_split), Intervals(c(-abs(sample_signal_split), abs(sample_signal_split)))))
      numerator_1 <- suppressWarnings(interval_intersection(interval_complement(phi_bounds_1), Intervals(c(-abs(sample_signal_1), abs(sample_signal_1)))))
      numerator_2<- suppressWarnings(interval_intersection(interval_complement(phi_bounds_2), Intervals(c(-abs(sample_signal_2), abs(sample_signal_2)))))



     

      ### SAVE ALL OF THE RESULTS.
      write(paste(c(beta,XORlev, cp, seed, "Tree-Values", "split", depth, 
                    p_split,
                    CI_split,
                    length(y1), length(y2), sample_signal_split,true_signal_split,NA,
                    sum(size(phi_bounds_split)), sum(size(numerator_split)),
                    full_rand_1, full_rand_2, full_rand_3), collapse=" "),
            file=filename,append=TRUE)
      
      write(paste(c(beta,XORlev, cp, seed, "Tree-Values", "child", depth, 
                    p_1,
                    CI_1,
                    length(y1), NA, sample_signal_1,true_signal_1,NA,
                    sum(size(phi_bounds_1)), sum(size(numerator_1)),
                    full_rand_1, full_rand_2, full_rand_3), collapse=" "),
            file=filename,append=TRUE)
      write(paste(c(beta,XORlev, cp, seed, "Tree-Values", "child", depth, 
                    p_2,
                    CI_2,
                    NA, length(y2), sample_signal_2,true_signal_2,NA,
                    sum(size(phi_bounds_2)), sum(size(numerator_2)),
                    full_rand_1, full_rand_2, full_rand_3), collapse=" "),
            file=filename,append=TRUE)
      j <- j+2

  }
  }
}

