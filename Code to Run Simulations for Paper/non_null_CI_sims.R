library(rpart)
library(intervals)
library(partykit)
library(rpart.utils)


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
  base_tree <- rpart::rpart(y~., data=dat,
    control=rpart.control(minbucket=minbucket,cp=cp, maxcompete=0,maxsurrogate=0, maxdepth=maxdepth), model=TRUE)


  ### If it build a non-trivial rpart tree,
  if (length(unique(base_tree$where)) > 1) {
  #### Go through EVERY SPLIT and EVERY NODE
  splitList <- getAllBranches(base_tree)
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
      CI_split <- computeCI(nu, y, sigma_y,phi_bounds_split, 0.05)
      CI_1 <- computeCI(nu1, y, sigma_y,phi_bounds_1, 0.05)
      CI_2 <- computeCI(nu2, y, sigma_y,phi_bounds_2, 0.05)

      ### Saving stuff about size of set S. Used to be used for plots.
      ### Didn't end up in any plots in final paper, but still saved bc interesting to look at.
      numerator_split <- suppressWarnings(interval_intersection(interval_complement(phi_bounds_split), Intervals(c(-abs(sample_signal_split), abs(sample_signal_split)))))
      numerator_1 <- suppressWarnings(interval_intersection(interval_complement(phi_bounds_1), Intervals(c(-abs(sample_signal_1), abs(sample_signal_1)))))
      numerator_2<- suppressWarnings(interval_intersection(interval_complement(phi_bounds_2), Intervals(c(-abs(sample_signal_2), abs(sample_signal_2)))))



      ### Naive Z INFERENCE. On the split and on the two individual children.
      seZ <- sqrt(sigma_y^2/length(y1)+sigma_y^2/length(y2))
      p_splitZ <- (1-pnorm(abs((sample_signal_split)/seZ)))*2
      p_1Z <- (1-pnorm(abs(mean(y1))/sqrt(sigma_y^2/length(y1))))*2
      p_2Z <- (1-pnorm(abs(mean(y2))/sqrt(sigma_y^2/length(y2))))*2
      CI_splitZ <- c(sample_signal_split - 1.96*seZ, sample_signal_split + 1.96*seZ)
      CI_1Z <- c(mean(y1) - 1.96*sigma_y/sqrt(length(y1)), mean(y1)  + 1.96*sigma_y/sqrt(length(y1)))
      CI_2Z <-  c(mean(y2) - 1.96*sigma_y/sqrt(length(y2)), mean(y2) + 1.96*sigma_y/sqrt(length(y2)))



      ### SAVE ALL OF THE RESULTS.
      write(paste(c(beta,XORlev, cp, seed, "Tree-Values", "split", depth, p_split, CI_split, length(y1), length(y2), sample_signal_split,true_signal_split,NA,
                    sum(size(phi_bounds_split)), sum(size(numerator_split))), collapse=" "),
            file=filename,append=TRUE)
      write(paste(c(beta,XORlev, cp, seed, "Tree-Values", "child", depth, p_1, CI_1, length(y1), NA, sample_signal_1,true_signal_1,NA,
                    sum(size(phi_bounds_1)), sum(size(numerator_1))), collapse=" "),
            file=filename,append=TRUE)
      write(paste(c(beta,XORlev, cp, seed, "Tree-Values", "child", depth, p_2, CI_2, NA, length(y2), sample_signal_2,true_signal_2,NA,
                    sum(size(phi_bounds_2)), sum(size(numerator_2))), collapse=" "),
            file=filename,append=TRUE)
      write(paste(c(beta,XORlev, cp, seed, "naiveZ", "split", depth, p_splitZ, CI_splitZ, length(y1), length(y2), sample_signal_split,true_signal_split,NA,
                   NA,NA), collapse=" "),
            file=filename,append=TRUE)
      write(paste(c(beta,XORlev, cp, seed, "naiveZ", "child", depth, p_1Z, CI_1Z, length(y1), NA, sample_signal_1,true_signal_1,NA,
                    NA,NA), collapse=" "),
            file=filename,append=TRUE)
      write(paste(c(beta,XORlev, cp, seed, "naiveZ", "child", depth, p_2Z, CI_2Z, NA, length(y2), sample_signal_2,true_signal_2,NA,
                    NA,NA), collapse=" "),
            file=filename,append=TRUE)


      j <- j+2

  }
  }

  ### Now repeat the process but use sample splitting!!!!!!!!!!!!
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
    j=1
    while (j < length(splitList)) {

      splits <- splitList[j][[1]]
      splitText1 <- paste(paste("dat$", splits), collapse=" & ")
      splitText2 <- paste(paste("dat$", splitList[j+1][[1]]), collapse=" & ")
      node1 <- eval(parse(text =splitText1))
      node2 <- eval(parse(text =splitText2))*2
      where <- node1+node2


      node1test <- node1[n2]
      node2test <- node2[n2]
      wheretest <- where[n2]
      mu_ytest <- mu_y[n2]

      y1 <- dat2[wheretest==1,]$y
      y2 <- dat2[wheretest==2,]$y
      mu1 <- mu_ytest[wheretest==1]
      mu2 <- mu_ytest[wheretest==2]

      nu <- (where==1)/sum(where==1) - (where==2)/sum(where==2)
      nu1 <- (where==1)/sum(where==1)
      nu2<- (where==2)/sum(where==2)




      true_sig_split <- as.numeric(t(nu)%*%mu_y)
      true_sig_1 <- as.numeric(t(nu1)%*%mu_y)
      true_sig_2 <- as.numeric(t(nu2)%*%mu_y)

      true_sig_split_test <- mean(mu1)-mean(mu2)
      true_sig_1_test <-  mean(mu1)
      true_sig_2_test <-  mean(mu2)

      samp_signal_split <-mean(y1)-mean(y2)
      samp_signal_1 <- mean(y1)
      samp_signal_2 <- mean(y2)

      depth = length(splits)
      se <- sqrt(sigma_y^2/length(y1)+sigma_y^2/length(y2))
      CI_split_split <- c(samp_signal_split - 1.96*se, samp_signal_split + 1.96*se)
      CI_1_split <- c(samp_signal_1 - 1.96*sqrt(sigma_y^2/length(y1)), samp_signal_1 + 1.96*sqrt(sigma_y^2/length(y1)))
      CI_2_split <- c(samp_signal_2 - 1.96*sqrt(sigma_y^2/length(y2)), samp_signal_2 + 1.96*sqrt(sigma_y^2/length(y2)))
      p_split_split <- (1-pnorm(abs((mean(y1)-mean(y2))/se)))*2
      p_1_split <- (1-pnorm(abs(mean(y1))/sqrt(sigma_y^2/length(y1))))*2
      p_2_split <- (1-pnorm(abs(mean(y2))/sqrt(sigma_y^2/length(y2))))*2


      j <- j+2

      write(paste(c(beta,XORlev, cp, seed, "splitting", "split", depth, p_split_split, CI_split_split, length(y1), length(y2), samp_signal_split,
                    true_sig_split, true_sig_split_test, NA,NA), collapse=" "),
            file=filename,append=TRUE)
      write(paste(c(beta,XORlev, cp, seed, "splitting", "child", depth, p_1_split, CI_1_split, length(y1), NA, samp_signal_1,
                    true_sig_1, true_sig_1_test,  NA, NA), collapse=" "),
            file=filename,append=TRUE)
      write(paste(c(beta,XORlev,cp, seed, "splitting", "child", depth, p_2_split, CI_2_split, NA, length(y2), samp_signal_2,
                    true_sig_2, true_sig_2_test,  NA,NA), collapse=" "),
            file=filename,append=TRUE)
    }
  }
}

