# ------------------------------------------------------------------------------#
# The following is an EM algorithm for probabilistic principal components       #
# analysis. Based on Tipping and Bishop, 1999, and also Murphy 2012             #
# Probabilistic ML, with some code snippets inspired by the ppca function used  #
# below. See also ModelFitting/EM Examples/EM for pca.R                         #
# ------------------------------------------------------------------------------#

#####################
### Main Function ###
#####################

PPCAEM = function(X, nComp=2, tol=.00001, maxits=100, showits=T){
  # Arguments X: numeric data, nComp: number of components
  # tol = tolerance level, maxits: maximum iterations, showits: show iterations
  require(pracma) # for orthonormal basis of W; pcaMethods package has also
  require(psych)  # for tr
  
  # starting points and other initializations
  Xorig = X
  X = X
  N = nrow(Xorig)
  D = ncol(Xorig)
  L = nComp
  NAs = is.na(Xorig)

  X[NAs] = 0
  S = (1/N) * t(X)%*%X
  evals = eigen(S)$values
  evecs = eigen(S)$vectors

  V = evecs[,1:L]
  Lambda = diag(evals[1:L])
  
  Z = t(replicate(L, rnorm(N)))                                    # latent variables
  sigma2 = 1/(D-L) * sum(evals[(L+1):D])                           # variance; average variance associated with discarded dimensions
  W = V %*% chol(Lambda-sigma2*diag(L)) %*% diag(L)                # loadings

  it = 0
  converged = FALSE
  ll = 0
  
  if (showits)                                                     # Show iterations
    cat(paste("Iterations of EM:", "\n"))
  while ((!converged) & (it < maxits)) {                  
    if(exists('W.new')){
      W.old = W.new
    }
    else {
      W.old = W
    }
    
    ll.old = ll
    
    proj = t(W.old%*%Z)
    Xnew = Xorig
    Xnew[NAs] = proj[NAs]
    X = Xnew
    
    Psi = sigma2*diag(L)
    M = t(W.old) %*% W.old + Psi
    
    W.new = S%*%W.old%*%solve(Psi + solve(M)%*%t(W.old)%*%S%*%W.old)   # E and M
    sigma2 = 1/D * tr(S - S%*%W.old%*%solve(M)%*%t(W.new))

    Z = solve(M)%*%t(W.new)%*%t(X)
    
  
    # log likelihood as in paper
#     ZZ = sigma2*solve(M) + Z%*%t(Z)
#     ll = .5*sigma2*D + .5*tr(ZZ) + .5*sigma2 * X%*%t(X) -
#          1/sigma2 * t(Z)%*%t(W.new)%*%t(X) + .5*sigma2 * tr(t(W.new)%*%W.new%*%ZZ)
#     ll = -sum(ll)
    
    # more straightforward
    ll = dnorm(X, mean=t(W.new%*%Z), sd=sqrt(sigma2), log=T)
    ll = -sum(ll)

    it = it + 1
    if (showits & (it == 1 | it%%5 == 0))                          # if showits, show first and every 5th iteration
      cat(paste(format(it), "...", "\n", sep = ""))
    converged = max(abs(ll.old-ll)) <= tol
  }
  
  W = orth(W.new)
  evs = eigen(cov(X %*% W))
  evecs = evs$vectors
  
  W = W %*% evecs
  Z = X %*% W
  Xrecon = Z %*% t(W)
  reconerr = sum((Xrecon-X)^2)
  
  if (showits)                                                     # Show last iteration
    cat(paste0(format(it), "...", "\n"))
  
  return(list(scores=Z, loadings=W, Xrecon=Xrecon, reconerr=reconerr, ll=ll, sigma2=sigma2))
}



###############
### Example ###
###############

### Set up data
# state.x77 is the data; various state demographics
X = scale(state.x77)
Xmiss = X

# create some missing values
set.seed(123)
NAindex = sample(length(X), 20)
Xmiss[NAindex] = NA

### run pca
outEM = PPCAEM(X=Xmiss, nComp=2, tol=1e-8, maxit=100)
outEM  

# Extract reconstructed values and loadings for comparison
Xrecon = outEM$Xrecon
loadingsEM = outEM$loadings
scoresEM = outEM$scores

# mean squared reconstruction error
mean((Xrecon-X)^2)  
mean((Xrecon[NAindex]-X[NAindex])^2)  


### compare to standard pca on full data set if desired
origpca =  princomp(scale(state.x77))
scores_origpca = origpca$scores[,1:2]
loadings_origpca = origpca$loadings[,1:2]
Xrecon_origpca = scores_origpca%*%t(loadings_origpca)

# examine difference from pca on complete data
# sum((abs(loadingsEM)-abs(origpca$loadings[,1:2]))^2)



#################################################
### compare results to output from pcaMethods ###
#################################################

### Run pca. Note that signs for loadings/scores may be opposite.
library(pcaMethods)
outpcam = pca(Xmiss, nPcs=2, threshold=1e-8, method='ppca', scale='none', center=F)
loadings_pcam = loadings(outpcam)
scores_pcam = scores(outpcam)

### compare loadings and scores
round(cbind(loadings_pcam, loadingsEM, loadings_origpca), 3)
sum((abs(loadings_pcam)-abs(loadingsEM))^2)
round(cbind(abs(scores_pcam), abs(scoresEM)), 2)

# compare reconstructed data sets
Xrecon_pcam = scores_pcam %*% t(loadings_pcam)
mean((Xrecon_pcam-X)^2)  
mean((Xrecon_pcam[NAindex]-X[NAindex])^2)  
mean(abs(Xrecon_pcam-Xrecon))

# plots
library(car)
scatterplotMatrix(cbind(X[,1], Xrecon[,1], Xrecon_pcam[,1]))
scatterplotMatrix(cbind(X[,2], Xrecon[,2], Xrecon_pcam[,2]))

ggplot2:::qplot(Xrecon[,1], Xrecon_pcam[,1], geom=c('point','smooth'))

scatterplotMatrix(cbind(scoresEM, scores_pcam), pch=19, cex=.75 )
