library(rpart)
library(tidyverse)
library(gridExtra)
library(patchwork)
library(treevalues)

#### ORIGINAL DATA
set.seed(35)
n <- 100
X1 <- rnorm(n, mean=0,sd=2)
X2 <- rnorm(n, mean=0,sd=2)
mu_y <- 10*(X1 < -1) - 10*(X1 > -1) + 15*(X1 > -1 & X2 > 2)
y <- rnorm(100, mean=mu_y, sd=1)

#### ORIGINAL TREE
tree1 <- rpart(y ~ X1+X2, control=rpart.control(maxdepth = 2, maxcompete=0, maxsurrogate=0, cp=0.0), model=TRUE)
nu <- (tree1$where==3)/sum((tree1$where==3)) - (tree1$where==4)/sum(tree1$where==4)
Pi_perp <- diag(rep(1,n)) - nu%*%t(nu)/sum(nu^2)
nureg <- (tree1$where==3)/sum((tree1$where==3))
Pi_perpreg <- diag(rep(1,n)) - nureg%*%t(nureg)/sum(nureg^2)


#### SIBLING NODE CASE: EXAMPLE PERTURBATIONS
splits <- getBranch(tree1,4)
library(intervals)
getInterval(tree1, nu, splits)

as.numeric(t(nu)%*%y)
sum(X1  <=  -0.89402)
sum(X2 < 1.988722)

#Panel 2- no change bc medium
phi2 <- 5
y2 <- as.numeric(Pi_perp%*%y + nu/sum(nu^2)*(phi2))
tree2 <- rpart(y2 ~ X1+X2, control=rpart.control(maxdepth = 2, maxcompete=0, maxsurrogate=0,cp=-1))

# Panel 3 - things change bc its 0 duh
phi0 <- 0
y0 <- as.numeric(Pi_perp%*%y + nu/sum(nu^2)*(phi0))
tree0 <- rpart(y0 ~ X1+X2, control=rpart.control(maxdepth = 2, maxcompete=0, maxsurrogate=0,cp=0))

# Panel 4- things change bc too big
phi4 <- 40
y4 <- as.numeric(Pi_perp%*%y + nu/sum(nu^2)*(phi4))
tree4 <- rpart(y4 ~ X1+X2, control=rpart.control(maxdepth = 2, maxcompete=0, maxsurrogate=0,cp=0))


#### SINGLE NODE CASE: EXAMPLE PERTURBATIONS
getInterval(tree1, nureg, splits)
getInterval(tree1, nureg, splits[2:1])

### THEY SWITCH ORDERS
phi4reg <-40
y4reg <- as.numeric(Pi_perpreg%*%y + nureg/sum(nureg^2)*(phi4reg))
tree4reg <- rpart(y4reg ~ X1+X2, control=rpart.control(maxdepth = 2, maxcompete=0, maxsurrogate=0))

### ACTUALLY ITS FINE AND IN DATASET
phi0reg <- 0
y0reg <- as.numeric(Pi_perpreg%*%y + nureg/sum(nureg^2)*(phi0reg))
tree0reg <- rpart(y0reg ~ X1+X2, control=rpart.control(maxdepth = 2, maxcompete=0, maxsurrogate=0, cp=-1))

##### 5 is the one not in dataset!!!!!
phi2reg <- 5
y2reg <- as.numeric(Pi_perpreg%*%y + nureg/sum(nureg^2)*(phi2reg))
tree2reg <- rpart(y2reg ~ X1+X2, control=rpart.control(maxdepth = 2, maxcompete=0, maxsurrogate=0, cp=-1))


###PLOTS FOR SIB

# Orig tree
t1 <- expression(paste("y = y'(", nu[sib]^T, "y,", nu[sib], ")"))
subtit1 <- expression(paste("=y'(-14.9,", nu[reg],")"))
t1full <- expression(paste("y = y'(", nu[sib]^T, "y,", nu[sib], ") = y'(-14.9,", nu[reg],")"))
p11 <- ggplot() + geom_point(data=NULL, aes(x=X1, y=X2, col=y)) +
  geom_vline(xintercept=tree1$splits[1,4], col="red")+
  geom_line(data=NULL, aes(x=seq(tree1$splits[1,4],5,length.out=1000), y=rep(tree1$splits[2,4],1000)),  col="red")+
  ylab("Feature 2")+
  geom_text(aes(x=5,y=3,label="R[B]"),parse=TRUE,cex=5)+
  geom_text(aes(x=5,y=1,label="R[A]"),parse=TRUE,cex=5)+
  ggtitle(t1full)+
  theme(axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=12),
        )


# Same as orig tree
t2 <- expression(paste("y'(5, ", nu[sib], ")"))
p12 <- ggplot() +
  geom_point(data=NULL, aes(x=X1, y=X2, col=y2)) +
  geom_vline(xintercept=tree2$splits[1,4], col="red")+
  geom_line(data=NULL, aes(x=c(tree2$splits[1,4],6), y=rep(tree2$splits[2,4],2)), col="red")+
  ggtitle(t2)+
  ylab("")
p12


# NOT same as orig tree but first split same
t3 <- expression(paste("y'(0, ", nu[sib], ")"))
p13 <- ggplot() + geom_point(data= NULL, aes(x=X1, y=X2, col=y0)) +
  geom_vline(xintercept=tree0$splits[1,4], col="red")+
  geom_line(data=NULL, aes(x=c(tree0$splits[1,4],6), y=rep(tree0$splits[2,4],2)), col="red")+
  ggtitle(t3)+
  guides(col=FALSE)+ylab("")+
  xlab("Feature 1")

# Not same as orig tree and splits change order!!!!!
t4 <- expression(paste("y'(40, ", nu[sib], ")"))
p14 <- ggplot() + geom_point(data=NULL, aes(x=X1, y=X2, col=y4)) +
  geom_hline(yintercept=tree4$splits[1,4], col="red")+
  geom_line(data=NULL, aes(y=c(tree4$splits[1,4],5.9), x=rep(tree4$splits[2,4],2)), col="red")+
  ggtitle(t4)+
  ylab("")+xlab("Feature 1")#+
  #scale_y_continuous(breaks=seq(-5, 5, 2.5))

### PLOTS FOR REGION

# Orig tree
t1reg <- expression(paste("y = y'(", nu[reg]^T, "y,", nu[reg], ")"))
subtit1reg <- expression(paste("=y'(-10,", nu[reg],")"))
t1fullreg <- expression(paste("y = y'(", nu[reg]^T, "y,", nu[reg], ") = y'(-10,", nu[reg],")"))
p21 <- ggplot() +
  geom_point(data=NULL, aes(x=X1, y=X2, col=y)) +
  geom_vline(xintercept=tree1$splits[1,4], col="red")+
  geom_line(data=NULL, aes(x=seq(tree1$splits[1,4],7,length.out=1000), y=rep(tree1$splits[2,4],1000)),  col="red")+
  ggtitle(t1fullreg)+
  ylab("Feature 2")+
  geom_text(aes(x=5,y=1,label="R[A]"),parse=TRUE,cex=5)

# NOT ORIG TREE. BOTH SPLITS ON X1
t22 <- expression(paste("y'(5, ", nu[reg], ")"))
p22 <- ggplot() + geom_point(data=NULL, aes(x=X1, y=X2, col=y2reg)) +
  geom_vline(xintercept=tree2reg$splits[1,4], col="red")+
  geom_vline(xintercept=tree2reg$splits[2,4], col="red")+
  ylab("")+
  ggtitle(t22)

t23 <- expression(paste("y'(0, ", nu[reg], ")"))
p23 <- ggplot() + geom_point(data=NULL, aes(x=X1, y=X2, col=y0reg)) +
  geom_vline(xintercept=tree0reg$splits[1,4], col="red")+
  geom_line(data=NULL, aes(x=c(tree0reg$splits[1,4],6), y=rep(tree0reg$splits[2,4],2)), col="red")+
  ggtitle(t23)+
  ylab("")


t24 <- expression(paste("y'(40, ", nu[reg], ")"))
p24 <- ggplot() +
  geom_point(data=NULL, aes(x=X1, y=X2, col=y4reg)) +
  geom_hline(yintercept=tree4reg$splits[1,4], col="red")+
  geom_line(data=NULL, aes(x=rep(tree4reg$splits[2,4],1000),
                           y=seq(-5, tree4reg$splits[1,4],length.out=1000)),
                            col="red")+
  ggtitle(t24)+
  ylab("")

#### MAKING THE PLOTS
colscheme <-scale_colour_gradientn(colours = hcl.colors(10, palette="viridis"),
                                   limits = range(y, y0, y4, y2, y4reg, y0reg),
                                   values = c(0, seq(0.3, 0.65, length.out = 8), 1),
                                   name = "y")
colscheme2 <-scale_colour_gradientn(colours = hcl.colors(12, palette="BuPu")[10:1],
                                   limits = range(y, y0, y4, y2, y4reg, y0reg),
                                   values = c(0, seq(0.3, 0.65, length.out = 8), 1),
                                   name = "y")

theme <-  theme_bw()+
  theme(panel.grid.major = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black"),
    panel.grid.minor = element_blank())+
  theme(axis.text.y=p11$theme$axis.text.y,
    axis.ticks.y = p11$theme$axis.ticks.y,
    axis.title.y=element_text(size=14),
    axis.title.x=element_text(size=14),
    axis.text.x=element_text(size=12),
    title=element_text(size=10))

p11 + p12 + p13 + p14 + p21 +p22+ p23+ p24 +
  plot_layout(guides="collect", nrow=2, byrow=TRUE) &
  xlim(-4,6) & ylim(-4,6) &
  #guides(col=FALSE) &
  coord_fixed() &
  theme &
  colscheme & xlab("Feature 1")

ggsave(
filename="~/Dropbox/Tree Values Paper/newformat_v30/Figures/squat_intuition.png")

