### The following code is based on algorithms noted in Murphy, 2012 Probabilistic Machine Learning.
### Specifically, Chapter 11, section 4.

###############################
### EM for gaussian mixture ###
###############################
gaussmixEM = function(params, X, clusters = 2, tol=.00001, maxits=100, showits=T){
  # Arguments are starting parameters (means, covariances, cluster probability), data, 
  # number of clusters desired, tolerance, maximum iterations, and whether to show iterations
  
  require(mvtnorm)
  # Starting points
  N = nrow(X)
  mu = params$mu
  var = params$var
  probs = params$probs
  
  # initializations
  ri = matrix(0, ncol=clusters, nrow=N)         # cluster 'responsibilities', i.e. probability of cluster membership for each observation i
  ll = 0                                        # loglike
  it = 0                                        # iteration count
  converged = FALSE                             # convergence
  
  if (showits)                                  # Show iterations if showits == true
    cat(paste("Iterations of EM:", "\n"))
  
  while (!converged & it < maxits) { 
    probsOld = probs
#   muOld = mu                                 # Use direct values or loglike for convergence
#   varOld = var
    llOld = ll
    riOld = ri
    
    ### E
    # Compute responsibilities
    for (k in 1:clusters){
      ri[,k] = probs[k] * dmvnorm(X, mu[k,], sigma = var[[k]], log=F)
    }
    ri = ri/rowSums(ri)
    
    ### M
    rk = colSums(ri)                            # rk is weighted average cluster membership size
    probs = rk/N
    for (k in 1:clusters){
      varmat = matrix(0, ncol=ncol(X), nrow=ncol(X))        # initialize to sum matrices
      for (i in 1:N){
        varmat = varmat + ri[i,k] * X[i,]%*%t(X[i,])
      }
      mu[k,] = (t(X) %*% ri[,k]) / rk[k]
      var[[k]] =  varmat/rk[k] - mu[k,]%*%t(mu[k,])
      ll[k] = -.5 * sum( ri[,k] * dmvnorm(X, mu[k,], sigma = var[[k]], log=T) )
    }
    ll = sum(ll)
    
    ### compare old to current for convergence
    parmlistold =  c(llOld, probsOld)           # c(muOld, unlist(varOld), probsOld)
    parmlistcurrent = c(ll, probs)              # c(mu, unlist(var), probs)
    it = it + 1
    if (showits & it == 1 | it%%5 == 0)         # if showits true, & it =1 or modulo of 5 print message
      cat(paste(format(it), "...", "\n", sep = ""))
    converged = min(abs(parmlistold - parmlistcurrent)) <= tol
  }
  
  clust = which(round(ri)==1, arr.ind=T)        # create cluster membership
  clust = clust[order(clust[,1]), 2]            # order accoring to row rather than cluster
  out = list(probs=probs, mu=mu, var=var, resp=ri, cluster=clust, ll=ll)
} 


#########################################
### Example 1: Old Faithful eruptions ###
#########################################

### This example uses Old Faithful geyser eruptions.  This is can be compared to the univariate code of EM Mixture.R
### See also http://www.geyserstudy.org/geyser.aspx?pGeyserNo=OLDFAITHFUL

# Create starting values
mustart = rbind(c(3, 60), c(3, 60.1))    # must be at least slightly different
covstart = list(cov(faithful), cov(faithful))
probs = c(.01, .99)

starts = list(mu=mustart, var=covstart, probs=probs)  # params is a list of mu var and probs 

# Run and examine
test = gaussmixEM(params=starts, X=as.matrix(faithful), clusters = 2, tol=1e-8, maxits=1500, showits=T)
str(test)

### graphical display
# library(ggplot2)
# qplot(x=eruptions, y=waiting, data=faithful) + ggtheme
# 
# ggplot(aes(x=eruptions, y=waiting), data=faithful) +
#   geom_point(aes(color=factor(test$cluster))) +
#   ggtheme
# 
# ggplot(aes(x=eruptions, y=waiting), data=faithful) +
#   geom_point(aes(color=test$resp[,1])) +
#   ggtheme
# 
# worst = apply(test$resp, 1, function(x) max(x) < .99) #relatively speaking, these are extremely well-separated clusters
# ggplot(aes(x=eruptions, y=waiting), data=faithful) +
#   geom_point(aes(color=worst)) +
#   ggtheme

### Compare to mclust results
library(mclust)
mclustmod = Mclust(faithful[,1:2], 2)
str(mclustmod,1)

# compare means
t(test$mu); mclustmod$parameters$mean

# compare variances
test$var; mclustmod$parameters$variance$sigma

# compare classifications, reverse in case arbitrary numbering of one of them is opposite
table(test$cluster, mclustmod$classification)
table(ifelse(test$cluster==2, 1, 2), mclustmod$classification)

# compare responsibilities; reverse one if arbitrary numbering of one of them is opposite
cbind(round(test$resp[,1], 2), round(mclustmod$z[,2], 2)) #cluster '1'
cbind(round(test$resp[,2], 2), round(mclustmod$z[,1], 2)) #cluster '2'

################################
### Example 2: Iris data set ###
################################
# Set up starting values
iris2 = iris[,1:4]

library(plyr)
mustart = daply(iris, 'Species', function(x) colMeans(x[,1:4])) + runif(4, 0, .5)    #add noise; function is notably sensitive to starts, but don't want to cheat too badly
covstart = dlply(iris, 'Species', function(x) var(x[,1:4]) + diag(runif(4, 0, .5)))
probs = c(.1, .2, .7)

starts = list(mu=mustart, var=covstart, probs=probs)  

# Run and examine
test = gaussmixEM(params=starts, X=as.matrix(iris2), clusters = 3, tol=1e-8, maxits=1500, showits=T)
table(test$cluster, iris$Species)

### Compare to mclust results
mclustIris = Mclust(iris[,1:4], 3)
table(mclustIris$classification, iris$Species)