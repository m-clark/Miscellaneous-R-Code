########################################################################################
### 'Noisy' gaussian process demo.  The matrix labeling is in keeping with Murphy    ###
### 2012 and Rasmussen and Williams 2006.  See those sources for more detail.        ###
### Murphy's matlab code can be found here: https://code.google.com/p/pmtk3/, though ###
### the relevant files are housed alongside this code.                               ###
###                                                                                  ###
### The goal of this code is to plot samples from the prior and posterior predictive ###
### of a gaussian process in which y = sin(x) + noise. It will reproduce an example  ###
### akin to figure 15.3 in Murphy 2012.                                              ###
########################################################################################


#################
### Functions ###
#################

# the mean function; in this case mean=0
muFn = function(x){
  x = sapply(x, function(x) x=0)
  x
}

# The covariance function; here it is the squared exponential kernel.
# l is the horizontal scale, sigmaf is the vertical scale, sigman the noise.
# See ?covSEiso in the gpr package for example, which is also based on Rasmussen and
# Williams Matlab code (gpml Matlab library)

Kfn = function(x, y=NULL, l=1, sigmaf=1, sigman=.5){
  if(!is.null(y)){
    sigmaf * exp( -(1/(2*l^2)) * as.matrix(dist(x, upper=T, diag=T)^2) ) + sigman*diag(length(x))    
  }
  else{
    sigmaf * exp( -(1/(2*l^2)) * as.matrix(dist(x, upper=T, diag=T)^2) )
  }  
}

#####################
### Preliminaries ###
#####################

l = 1           # for l, sigmaf, sigman, see note at covariance function
sigmaf = 1      
sigman = .25 
keps = 1e-8     # see note at Kstarstar
nprior = 5      # number of prior draws
npostpred = 3   # number of posterior predictive draws

##################
### Prior plot ###
##################

### data setup
require(MASS)
xg1 = seq(-5, 5, .2)
yg1 = mvrnorm(nprior, mu=muFn(xg1), Sigma=Kfn(xg1, l=l, sigmaf=sigmaf, sigman=sigman)) 

### plot
library(ggplot2); library(reshape2)

# reshape data for plotting
gdat = melt(data.frame(x=xg1, y=t(yg1), sd=apply(yg1, 2, sd)), id=c('x', 'sd'))
# head(gdat) # inspect if desired

g1 = ggplot(aes(x=x, y=value), data=gdat) + 
  geom_line(aes(group=variable), color='#FF5500', alpha=.5) +
  labs(title='Prior') +
  ggtheme

# g1 # show plot

####################################
### generate noisy training data ###
####################################

Xtrain = 15*(runif(20)-.5)  
nTrain = length(Xtrain)
ytrain = sin(Xtrain) + rnorm(n=nTrain, sd=.1)  # kept sine function for comparison to noise free result

Xtest = seq(-7.5, 7.5, length=200)
nTest = length(Xtest)

#####################################
### generate posterior predictive ###
#####################################

### Create Ky, K*, and K** matrices as defined in the texts
Ky = Kfn(x=Xtrain, y=ytrain, l=l, sigmaf=sigmaf, sigman=sigman)
K_ = Kfn(c(Xtrain, Xtest), l=l, sigmaf=sigmaf, sigman=sigman)                    # initial matrix
Kstar = K_[1:nTrain, (nTrain+1):ncol(K_)]                                        # dim = N x N*
tKstar = t(Kstar)                                                                # dim = N* x N
Kstarstar = K_[(nTrain+1):nrow(K_), (nTrain+1):ncol(K_)] + keps*diag(nTest)      # dim = N* x N*; the keps part is for positive definiteness
Kyinv = solve(Ky)

# calculate posterior mean and covariance
postMu = muFn(Xtest) + tKstar %*% Kyinv %*% (ytrain-muFn(Xtrain))
postCov = Kstarstar - tKstar %*% Kyinv %*% Kstar
s2 = diag(postCov)
# R = chol(postCov)  
# L = t(R)      # L is used in alternative formulation below based on gaussSample.m

# generate draws from posterior predictive
y2 = data.frame(t(mvrnorm(npostpred, mu=postMu, Sigma=postCov)))
# y2 = data.frame(replicate(npostpred, postMu + L %*% rnorm(postMu))) # alternative

#################################
### Posterior predictive plot ###
#################################

# reshape data for plotting
gdat = melt(data.frame(x=Xtest, y=y2, fmean=postMu, selower=postMu-2*sqrt(s2), seupper=postMu+2*sqrt(s2)),
            id=c('x', 'fmean','selower', 'seupper'))

g2 = ggplot(aes(x=x, y=value), data=gdat) + 
  geom_ribbon(aes(ymin=selower, ymax=seupper,group=variable), fill='gray90') +
  geom_line(aes(group=variable), color='#FF5500', alpha=.5) +
  geom_line(aes(group=variable, y=fmean), color='navy') +
  geom_point(aes(x=Xtrain, y=ytrain), data=data.frame(Xtrain, ytrain)) +
  labs(title='Posterior Predictive') +
  ggtheme

# g2

####################################################
### Plot prior and posterior predictive together ###
####################################################

library(gridExtra)
grid.arrange(g1, g2, ncol=2)