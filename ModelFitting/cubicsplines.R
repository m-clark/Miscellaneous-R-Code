#################################################
### Cubic Spline example. See Wood (2006)     ###
### Generalized Additive Models or my handout ###
### at  www.nd.edu/~mclark19/learn/GAMS.pdf   ###
#################################################
 
### Create the data
size = c(1.42,1.58,1.78,1.99,1.99,1.99,2.13,2.13,2.13,
         2.32,2.32,2.32,2.32,2.32,2.43,2.43,2.78,2.98,2.98)
wear = c(4.0,4.2,2.5,2.6,2.8,2.4,3.2,2.4,2.6,4.8,2.9,
         3.8,3.0,2.7,3.1,3.3,3.0,2.8,1.7)
x = size-min(size)
x = x/max(x)
d = data.frame(wear, x)

### cubic spline function
rk <- function(x, z) {
  ((z-0.5)^2 - 1/12) * ((x-0.5)^2 - 1/12)/4 -
    ((abs(x-z)-0.5)^4 - (abs(x-z)-0.5)^2/2 + 7/240) / 24
}

### Generate the model matrix
splX <- function(x, knots){
  q <- length(knots) + 2 # number of parameters
  n <- length(x) # number of observations
  X <- matrix(1, n, q) # initialized model matrix
  X[,2] <- x # set second column to x
  X[,3:q] <- outer(x, knots, FUN=rk) # remaining to cubic spline basis
  X
}

splS <- function(knots) {
  q = length(knots) + 2
  S = matrix(0, q, q) # initialize matrix
  S[3:q,3:q] = outer(knots, knots, FUN=rk) # fill in non-zero part
  S
}

### Matrix square root function; Note that there are various packages with their own.
matSqrt <- function(S){
  d = eigen(S, symmetric=T)
  rS = d$vectors %*% diag(d$values^.5) %*% t(d$vectors)
  rS
}

### Penalized fitting function
prsFit <- function(y, x, knots, lambda){
  q = length(knots) + 2 # dimension of basis
  n = length(x) # number of observations
  Xa = rbind(splX(x, knots), matSqrt(splS(knots))*sqrt(lambda)) # augmented model matrix
  y[(n+1):(n+q)] = 0 #augment the data vector
  lm(y ~ Xa-1) # fit and return penalized regression spline
}


#################
### Example 1 ###
#################

# Unpenalized
knots = 1:4/5
X = splX(x, knots) # generate model matrix
mod1 = lm(wear~X-1) # fit model
xp <- 0:100/100 # x values for prediction
Xp <- splX(xp, knots) # prediction matrix

# Base R plot
plot(x, wear, xlab='Scaled Engine size', ylab='Wear Index', pch=19,
     col="#FF8000", cex=.75, col.axis='gray50', bty='n')
lines(xp, Xp%*%coef(mod1), col='#2957FF') 

# ggplot
# library(ggplot2)
# ggplot(aes(x=x, y=wear), data=data.frame(x,wear))+
#   geom_point(color="#FF8000") +
#   geom_line(aes(x=xp, y=Xp%*%coef(mod1)), data=data.frame(xp,Xp), color="#2957FF") +
#   ggtheme

#################
### Example 2 ###
#################

# Add penalty lambda
knots = 1:7/8
d2 = data.frame(x=xp)
for (i in c(.1, .01, .001, .0001, .00001, .000001)){
  mod2 = prsFit(y=wear, x=x, knots=knots, lambda=i) #fit penalized regression
  #spline choosing lambda
  Xp = splX(xp, knots) #matrix to map parameters to fitted values at xp
  d2[,paste0('lambda = ', i)] = Xp%*%coef(mod2)
}

# Examine
# head(d2)

### Visualize via base R plot
par(mfrow=c(2,3))
for (i in c(.1, .01, .001, .0001, .00001, .000001)){
  mod2 = prsFit(y=wear, x=x, knots=knots, lambda=i)
  Xp = splX(xp, knots)
  plot(x, wear, xlab='Scaled Engine size', main=paste0('lambda = ', i), pch=19,
       col="#FF8000", cex=.75, col.axis='gray50', bty='n')
  lines(xp, Xp%*%coef(mod2), col='#2957FF')
}

### Visualize via ggplot
# library(ggplot2); library(reshape2)
# d3 = melt(d2, id='x')
# ggplot(aes(x=x, y=wear), data=d) +
#   geom_point(col='#FF8000') +
#   geom_line(aes(x=x,y=value), col="#2957FF", data=d3) +
#   facet_wrap(~variable) +
#   ggtheme

