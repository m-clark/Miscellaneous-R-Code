### The following code is based on algorithms noted in Murphy, 2012 Probabilistic Machine Learning.
### Specifically, Chapter 11, section 4.

###############################
### EM for gaussian mixture ###
###############################
gaussmixEM = function(params, X, clusters = 2, tol=.00001, maxits=100, showits=T){     
  # Arguments are starting parameters (means, covariances, cluster probability), data, number of clusters desired, tolerance,
  # maximum iterations, and whether to show iterations
  
  # Starting points
  N = nrow(X)
  nams = names(params)
  mu = params$mu
  var = params$var
  probs = params$probs
  
  # Other initializations
  ri = matrix(0, ncol=clusters, nrow=N)  #initialize cluster 'responsibilities', i.e. probability of cluster membership for each observation i
  it = 0
  converged = FALSE
  
  if (showits)                                  # Show iterations
    cat(paste("Iterations of EM:", "\n"))
  
  while ((!converged) & (it < maxits)) { 
    probsOld = probs
    muOld = mu
    varOld = var
    riOld = ri
    
    ### E
    # Compute responsibilities
    for (k in 1:clusters){
      ri[,k] = probs[k] * dnorm(X, mu[k], sd = sqrt(var[k]), log=F)
    }
    ri = ri/rowSums(ri)
    
    ### M
    rk = colSums(ri)                           # rk is the weighted average cluster membership size
    probs = rk/N
    mu = (t(X) %*% ri) / rk                      # could do mu and var via log likelihood here but this is more straightforward
    var = (t(X^2) %*% ri) / rk - mu^2

    parmlistold = rbind(probsOld, muOld, varOld)
    parmlistcurrent = rbind(probs, mu, var)
    it = it + 1
    if (showits & it == 1 | it%%5 == 0)        # if showits true, & it =1 or divisible by 5 print message
      cat(paste(format(it), "...", "\n", sep = ""))
    converged = max(abs(parmlistold - parmlistcurrent)) <= tol
  }
  
  clust = which(round(ri)==1, arr.ind=T)       # create cluster membership
  clust = clust[order(clust[,1]), 2]           # order accoring to row rather than cluster
  out = list(probs=probs, mu=mu, var=var, resp=ri, cluster=clust)
} 


### This example uses Old Faithful geyser eruptions.  This is only a univariate mixture for either eruption time or wait time.
### The next will be doing both variables, i.e. multivariate normal.  'Geyser' is supposedly more accurate, though seems to have 
### arbitrarily assigned some duration values.  See also http://www.geyserstudy.org/geyser.aspx?pGeyserNo=OLDFAITHFUL, but that only has
### intervals. Some July 1995 data is available

### faithful data set
data(faithful)
head(faithful)

# starting parameters; requires mean, variance and class probabilitiy
params1 = list(mu=c(2, 5), var=c(1, 1), probs=c(.5, .5))  # note that starts from mean must be in data range or it will break.  
params2 = list(mu=c(50, 90), var=c(1, 15), probs=c(.5, .5))  

X1 = matrix(faithful[,1])
X2 = matrix(faithful[,2])

test1 = gaussmixEM(params1, X=X1, tol=1e-8)
test2 = gaussmixEM(params2, X=X2, tol=1e-8)

# Compare to flexmix package results
library(flexmix)
flexmod1 = flexmix(X1~1, k=2, control=list(tolerance=1e-8, iter.max=100))
flexmod2 = flexmix(X2~1, k=2, control=list(tolerance=1e-8, iter.max=100))

# the following provides means, variances and probability of group membership;
# note that the cluster is arbitrary so cluster 1 for one model may be cluster 2 in another

### Eruptions
meanvar = rbind(test1$mu, sqrt(test1$var)); rownames(meanvar) = c('means', 'variances')
meanvarFlex = parameters(flexmod1); rownames(meanvarFlex) = c('means', 'variances')
probMembership = test1$probs
probMembershipFlex = flexmod1@size/sum(flexmod1@size)

list(params=cbind(meanvarFlex, meanvar), clusterpobs=cbind(probMembership, probMembershipFlex) )

### Waiting
meanvar = rbind(test2$mu, sqrt(test2$var)); rownames(meanvar) = c('means', 'variances')
meanvarFlex = parameters(flexmod2); rownames(meanvarFlex) = c('means', 'variances')
probMembership = test2$probs
probMembershipFlex = flexmod2@size/sum(flexmod2@size)

list(params=cbind(meanvarFlex, meanvar), clusterpobs=cbind(probMembership, probMembershipFlex) )

# MASS version (reversed columns)
# These don't look even remotely the same data on initial inspection; 'geyser' is even more rounded and of opposite conclusion;
# Turns out geyser is offset by 1, such that duration 1 should be coupled with waiting 2 and on down
# Still the rounding at 2 and 4 (and whatever division was done on duration) makes this fairly poor data

library(MASS)
geyser = data.frame(duration=geyser$duration[-299], waiting=geyser$waiting[-1])

# compare to faithful
layout(1:2); plot(faithful); plot(geyser)

X3 = matrix(geyser[,1]) 
X4 = matrix(geyser[,2])


### MASS version
test3 = gaussmixEM(params1, X=X3, tol=1e-8)
test4 = gaussmixEM(params2, X=X4, tol=1e-8)

flexmod3 = flexmix(X3~1, k=2, control=list(tolerance=1e-8, iter.max=100))
flexmod4 = flexmix(X4~1, k=2, control=list(tolerance=1e-8, iter.max=100))

# note variability differences compared to faithful dataset
# Eruptions/Duration
meanvar = rbind(test3$mu, sqrt(test3$var)); rownames(meanvar) = c('means', 'variances')
meanvarFlex = parameters(flexmod3); rownames(meanvarFlex) = c('means', 'variances')
probMembership = test3$probs
probMembershipFlex = flexmod3@size/sum(flexmod3@size)

list(params=cbind(meanvarFlex, meanvar), clusterpobs=cbind(probMembership, probMembershipFlex) )

# Waiting
meanvar = rbind(test4$mu, sqrt(test4$var)); rownames(meanvar) = c('means', 'variances')
meanvarFlex = parameters(flexmod4); rownames(meanvarFlex) = c('means', 'variances')
probMembership = test4$probs
probMembershipFlex = flexmod4@size/sum(flexmod4@size)

list(params=cbind(meanvarFlex, meanvar), clusterpobs=cbind(probMembership, probMembershipFlex) )

### Some plots; ggtheme available at https://github.com/m-clark/Miscellaneous-R-Code/blob/master/Other/ggtheme.R
library(ggplot2)
qplot(x=eruptions, y=waiting, data=faithful) + ggtheme

ggplot(aes(x=eruptions, y=waiting), data=faithful) +
  geom_point(aes(color=factor(test2$cluster))) +
  ggtheme

ggplot(aes(x=eruptions, y=waiting), data=faithful) +
  geom_point(aes(color=test2$resp[,1])) +
  ggtheme
