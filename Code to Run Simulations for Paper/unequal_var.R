library(rpart)
library(treevalues)


setwd("~/Dropbox/Tree Values Paper/Code Reviewer Response")
n <- 200
p <-10
sigma_y_true <- c(rep(5, n/2), rep(10, n/2))
nTrials <- 500
pvals <- rep(0, nTrials)
CIs <- matrix(0, nrow=nTrials, ncol=2)
truths <- rep(0, nTrials)
filename<- "unequal_var.csv"

for (i in 1:nTrials) {
  set.seed(i)
  print(i)

  X <- MASS::mvrnorm(n, rep(0,p), diag(rep(1,p)))
  y <- rnorm(n, 0, sigma_y_true)
  #y <- rpois(n,10)
  sigma_y <- sd(y)
  dat <- data.frame(y=y,X=X)
  nameX <- sapply(1:p, function(u) paste0("X",u))
  names(dat) = c("y", nameX)

  ### How did I pick this initially? Lol. 
  cp = 50/sum((y-mean(y))^2)

  base_tree <- rpart::rpart(y~., data=dat,
    control=rpart.control(minbucket=1,cp=cp, maxcompete=0,maxsurrogate=0, maxdepth=3), model=TRUE)

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
      nu <- (where==1)/sum(where==1) - (where==2)/sum(where==2)
    
      sample_signal_split <- as.numeric(t(nu)%*%y)
   
      phi_bounds_split <- getInterval(base_tree, nu,splits)
     
      depth = length(splits)
      p_split <- correctPVal(phi_bounds_split, nu, y, sigma_y)
      #p_1 <- correctPVal(phi_bounds_1, nu1, y, sigma_y)
      #p_2 <- correctPVal(phi_bounds_2, nu2, y, sigma_y)
      se <- sqrt(sigma_y^2/length(y1)+sigma_y^2/length(y2))
      p_split_naive <- 2*(1-pnorm(abs(sample_signal_split), mean=0, sd=se))
      #p_1_naive <-  2*(1-pnorm(abs(sample_signal_1), mean=0, sd=sqrt(sigma_y^2/length(y1))))
      #p_2_naive <-  2*(1-pnorm(abs(sample_signal_2), mean=0, sd=sqrt(sigma_y^2/length(y2))))
      j <- j+2

      write(paste(c(cp, i, "Tree-Values", "split", depth, p_split, length(y1), length(y2), sample_signal_split,'\n'), collapse=" "),
            file=filename,append=TRUE)
      
      write(paste(c(cp, i, "NaiveZ", "split", depth, p_split_naive, length(y1), length(y2), sample_signal_split, '\n'), collapse=" "),
            file=filename,append=TRUE)
  }
  } 

  ###### SAMPLE SPLITTING!!!!!!!!!!!!!!!!
  n1 <- sample(1:n, size=floor(n/2))
  n2 <- setdiff(1:n, n1)
  dat1 <- dat[n1,]
  dat2 <- dat[-n1,]
  cp = 50/sum((dat1$y-mean(dat1$y))^2)

  split_tree <- rpart::rpart(y~., data=dat1, control=rpart.control(minbucket=1,cp=cp, maxcompete=0,maxsurrogate=0, maxdepth=3), model=TRUE)

  split_tree$frame$yval = 1:NROW(split_tree$frame)
  predict = predict(split_tree, newdata=dat2)
  full_predict = predict(split_tree, newdata=dat)
  test_predict = predict(split_tree, newdata=dat2)
  terminalNodes <- sort(unique(full_predict))

  if (length(terminalNodes) > 1) {
    splitList <- getBranch(split_tree)
    j=1

    ### DEFINITIONS OF DETECTION ARE ON ALL DATA: train + test
    while (j < length(splitList)) {
      splits <- splitList[j][[1]]
      splitText1 <- paste(paste("dat$", splits), collapse=" & ")
      splitText2 <- paste(paste("dat$", splitList[j+1][[1]]), collapse=" & ")
      j <- j+2
      node1 <- eval(parse(text =splitText1))
      node2 <- eval(parse(text =splitText2))*2
      where <- node1+node2
      where_test <- (node1+node2)[n2]

      ##### DO INFERENCE ON TEST SET
      y1 <- dat2[where_test==1,]$y
      y2 <- dat2[where_test==2,]$y

      nu_full <- (where==1)/sum(where==1) - (where==2)/sum(where==2)
      sample_signal_split <- mean(y1)-mean(y2)
      se <- sqrt(sigma_y^2/length(y1) + sigma_y^2/length(y2))
      p_split <- (1-pnorm(abs((mean(y1)-mean(y2))/se)))*2
      depth = length(splits)

      write(paste(c(cp, i, "splitting", "split", depth, p_split, length(y1), length(y2), sample_signal_split, '\n'), collapse=" "),
            file=filename,append=TRUE)
    }
  }
}

