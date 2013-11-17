### 'noise free' gaussian process demo.  The matrix labeling is in keeping with Murphy  ###
### 2012 and Rasmussen and Williams 2006.  See those sources for more detail.  Murphy's ###
### matlab code can be found here: https://code.google.com/p/pmtk3/; specifically it is ###
### the file gprDemoNoiseFree.  Note that code refers to a sq_dist function that is not ###
### available in Octave.


### The goal of this code is to plot samples from the prior and posterior predictive 
### of a gaussian process in which y = sin(x)


###################################
### Initial setup and functions ###
###################################

l = 1
sigmaf = 1
keps = 1e-8  

# the mean function; in this case mean=0
muFn = function(x){
  x = sapply(x, function(x) x=0)
  x
}

# the covariance function; here it is the squared exponential kernel
# l is the horizontal scale, sigmaf is the vertical scale.
# see ?covSEiso in the gpr package for example, which is also based on Rasmussen and
# Williams Matlab code (gpml library); includes the missing sq_dist function

Kfn = function(x, l=1, sigmaf=1){
    sigmaf * exp( -(1/(2*l^2)) * as.matrix(dist(x, upper=T, diag=T)^2) )
}

# data setup
require(MASS)
x = seq(-5, 5, .2)
y = mvrnorm(5, mu=muFn(x), Sigma=Kfn(x, l=l, sigmaf=sigmaf)) 

# plot prior
library(ggplot2); library(reshape2)
# reshape data for plotting
gdat = melt(data.frame(x, y=t(y), sd=apply(y, 2, sd)), id=c('x', 'sd'))
# head(gdat) # inspect if desired

g1 = ggplot(aes(x=x, y=value), data=gdat) + 
  geom_line(aes(group=variable), color='#FF5500') +
  labs(title='Prior') +
  ggtheme

# g1

#########################################
### generate noise-less training data ###
#########################################
Xtrain = c(-4, -3, -2, -1, 1)
ftrain = sin(Xtrain)
nTrain = length(Xtrain)
nTest = length(x)

#####################################
### generate posterior predictive ###
#####################################
# Create K, K*, and K** matrices as defined in the texts
K = Kfn(Xtrain)  
Kstar = Kfn(c(Xtrain,x))[-c(1:nTrain), 1:nTrain]
Kstar = t(Kstar)
Kstarstar = Kfn(x) + keps*diag(nTest)
Kinv = solve(K)

# calculate posterior mean and covariance
postMu = muFn(x) + t(Kstar) %*% Kinv %*% (ftrain-muFn(Xtrain))
postCov = Kstarstar - t(Kstar) %*% Kinv %*% Kstar
s2 = diag(postCov)
R = chol(postCov)
L = t(R)

# generate draws from posterior predictive
# y2 = data.frame(replicate(3, postMu + L %*% rnorm(postMu))) # alternative
y2 = data.frame(t(mvrnorm(3, mu=postMu, Sigma=postCov)))

### Plot
gdat = melt(data.frame(x, y=y2, selower=postMu-2*sqrt(s2), seupper=postMu+2*sqrt(s2)),
            id=c('x', 'selower', 'seupper'))

g2 = ggplot(aes(x=x, y=value), data=gdat) + 
  geom_ribbon(aes(ymin=selower, ymax=seupper,group=variable), fill='gray90') +
  geom_line(aes(group=variable), color='#FF5500') +
  geom_point(aes(x=Xtrain, y=ftrain), data=data.frame(Xtrain, ftrain)) +
  labs(title='Posterior Predictive') +
  ggtheme

# g2

####################################################
### Plot prior and posterior predictive together ###
####################################################
library(gridExtra)
grid.arrange(g1, g2, ncol=2)