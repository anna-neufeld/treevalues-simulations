library(rpart)
library(treevalues)


setwd("~/Dropbox/Tree Values Paper/Code Reviewer Response/")
n <- 200
p <-10
sigma_y <- 5
nTrials <- 1000
filename<- "null_res_var_est.csv"

for (i in 1:nTrials) {
  set.seed(i)
  print(i)

  X <- MASS::mvrnorm(n, rep(0,p), diag(rep(1,p)))
  y <- rnorm(n, 0, sigma_y)
  dat <- data.frame(y=y,X=X)
  nameX <- sapply(1:p, function(u) paste0("X",u))
  names(dat) = c("y", nameX)

  cp = 200/sum((y-mean(y))^2)

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
      p_split_cons <- correctPVal(phi_bounds_split, nu, y, sd(y))
      p_split_SSE <- correctPVal(phi_bounds_split, nu, y, sqrt(sum((y-predict(base_tree))^2/(n-length(unique(base_tree$where))))))
      p_split_med_yc <- correctPVal(phi_bounds_split, nu, y, sqrt(median((y - median(y))^2)/qchisq(0.5, df=1)))
      p_split_mad <- correctPVal(phi_bounds_split, nu, y, median(abs(y - median(y))))
      
      
      j <- j+2

      write(paste(c(cp, i, "known", "split", depth, p_split, length(y1), length(y2), sample_signal_split,'\n'), collapse=" "),
            file=filename,append=TRUE)
      write(paste(c(cp, i, "cons", "split", depth, p_split_cons, length(y1), length(y2), sample_signal_split,'\n'), collapse=" "),
            file=filename,append=TRUE)
      write(paste(c(cp, i, "SSE", "split", depth, p_split_SSE, length(y1), length(y2), sample_signal_split,'\n'), collapse=" "),
            file=filename,append=TRUE)
      write(paste(c(cp, i, "med_yc", "split", depth, p_split_med_yc, length(y1), length(y2), sample_signal_split,'\n'), collapse=" "),
            file=filename,append=TRUE)
      write(paste(c(cp, i, "mad", "split", depth, p_split_mad, length(y1), length(y2), sample_signal_split,'\n'), collapse=" "),
            file=filename,append=TRUE)
    
      
  }
  } 
}

