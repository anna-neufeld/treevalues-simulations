library(rpart)
library(tidyverse)
library(rpart.plot)


### Generate data, build initial tree
set.seed(18)
n=100
X <- MASS::mvrnorm(n, rep(0,2), diag(rep(10,2)))
y <- -rnorm(n, 0, 1)
dat <- data.frame(y=y,X=X)
nameX <- sapply(1:2, function(u) paste0("X",u))
names(dat) = c("y", nameX)
base_tree <- rpart(y~X, model=TRUE,
                          control=rpart.control(maxdepth=1,
                                                maxcompete=0,
                                               maxsurrogate=0 ))


### Original Scatter Plot, all data
ggplot(data=NULL, aes(x=X[,1], y=X[,2], col=y)) + geom_point()+
           xlab(expression(X[1])) + ylab(expression(X[2])) + theme_bw() +
             theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+
             theme(axis.title=element_text(size=20), legend.text = element_text(size=16),
              legend.title=element_text(size=20),
             axis.text=element_text(size=20))+ coord_fixed()
ggsave("~/panel1.png")

### rpart plot
rpart.plot(base_tree, extra=1)



### Original Scatter Plot, draw tree split
ggplot(data=NULL, aes(x=X[,1], y=X[,2], col=y)) + geom_point()+
  xlab(expression(X[1])) + ylab(expression(X[2])) + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+
  geom_vline(xintercept = base_tree$splits[4], col="red")+
  theme(axis.title=element_text(size=20), legend.text = element_text(size=16),
        legend.title=element_text(size=20),
        axis.text=element_text(size=20))+ coord_fixed()
ggsave("~/panel3.png")



### pvalue
sigma_y <- 1
y1 <- y[base_tree$where==2]
y2 <- y[base_tree$where==3]
se <- sqrt(sigma_y^2/length(y1) + sigma_y^2/length(y2))
pnaive <- (1-pnorm(abs((mean(y1)-mean(y2))/se)))*2
print(pnaive)





### Sample splitting
set.seed(2)
train <- sample(1:n, size=n/2)
test <- (1:n)[-train]
split_tree <- rpart(y~., data=dat[train,], model=TRUE,
                    control=rpart.control(maxdepth=1,
                                          maxcompete=0,
                                          maxsurrogate=0))

### rpart plot of sample splitting tree
rpart.plot(split_tree, extra=1)

### Scatter plot of training data
ggplot(data=NULL, aes(x=X[train,1], y=X[train,2], col=y[train])) + geom_point()+
  xlab("Feature 1") + ylab("Feature 2") + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+
  theme(axis.title=element_text(size=22), legend.text = element_text(size=16),
        legend.title=element_text(size=14),
        axis.text=element_text(size=20))+ coord_fixed()+
  labs(col="y")+xlim(-8,8)+
  ylim(-8,8)

### Scatter plot with tree line added
ggplot(data=NULL, aes(x=X[train,1], y=X[train,2], col=y[train])) + geom_point()+
  xlab("Feature 1") + ylab("Feature 2") + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+
  geom_hline(yintercept = split_tree$splits[4], col="red")+
  theme(axis.title=element_text(size=22), legend.text = element_text(size=16),
        legend.title=element_text(size=14),
        axis.text=element_text(size=20))+ coord_fixed()+labs(col="y")+xlim(-8,8)+
  ylim(-8,8)

### Scatter plot of TEST SET with tree line added
ggplot(data=NULL, aes(x=X[test,1], y=X[test,2], col=y[test])) + geom_point()+
  xlab("Feature 1") + ylab("Feature 2") + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+
  geom_hline(yintercept = split_tree$splits[4], col="red")+
  theme(axis.title=element_text(size=22), legend.text = element_text(size=16),
        legend.title=element_text(size=14),
        axis.text=element_text(size=20))+ coord_fixed()+labs(col="y")+xlim(-8,8)+
  ylim(-8,8)

### pvalue
sigma_y <- 1

## (getting predictions on test set)
temp_tree <- split_tree
temp_tree$frame$yval <- rownames(split_tree$frame)
testpreds <- predict(temp_tree, newdata=dat[test,])
y1 <- y[test][testpreds==2]
y2 <- y[test][testpreds==3]
se <- sqrt(sigma_y^2/length(y1) + sigma_y^2/length(y2))
pnaive <- (1-pnorm(abs((mean(y1)-mean(y2))/se)))*2
print(pnaive)


## Note that this code literally comes from their VISTREE github page.
covariates <-blsdata %>% dplyr::select(hunger, disinhibition, resteating, rrvfood, liking, wanting, kcal24h0)
covariates<-data.frame(scale(covariates))
p4<-rpart(kcal24h0~., model = TRUE, data = covariates,control = rpart.control(cp = 0.0185, maxcompete = 0, maxsurrogate = 0))
treeval.plot(p4, digits=3,nn=FALSE,printn=FALSE,inferenceType=1, space=1.7)





