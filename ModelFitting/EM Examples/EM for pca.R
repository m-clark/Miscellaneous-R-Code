# --------------------------------------------------------------------#
# The following is an EM algorithm for principal components analysis. #
# See Murphy, 2012 Probabilistic Machine Learning 12.2.5. Some of the #
# constructed object is based on output from pca function used below. # 
# --------------------------------------------------------------------#

#####################
### Main Function ###
#####################

PCAEM = function(X, nComp=2, tol=.00001, maxits=100, showits=T){
  # Arguments X: numeric data, nComp: number of components
  # tol = tolerance level, maxits: maximum iterations, showits: show iterations
  require(pracma) # for orthonormal basis of W; pcaMethods package has also
  
  # starting points and other initializations
  N = nrow(X)
  D = ncol(X)
  L = nComp
  Xt = t(X)
  Z = t(replicate(L, rnorm(N)))                                    # latent variables
  W = replicate(L, rnorm(D))                                       # loadings
  
  it = 0
  converged = FALSE
    
  if (showits)                                                     
    cat(paste("Iterations of EM:", "\n"))
  while ((!converged) & (it < maxits)) {                           # while no convergence and we haven't reached our max iterations do this stuff
    Z.old = Z                                                      # create 'old' values for comparison
    Z = solve(t(W)%*%W) %*% crossprod(W, Xt)                       # E
    W = Xt%*%t(Z) %*% solve(tcrossprod(Z))                         # M

    it = it + 1
    if (showits & (it == 1 | it%%5 == 0))                          # if showits, show first and every 5th iteration
      cat(paste(format(it), "...", "\n", sep = ""))
    converged = max(abs(Z.old-Z)) <= tol
  }
  
  # calculate reconstruction error
  Xrecon = W %*% Z
  reconerr = sum((Xrecon-t(X))^2)
  
  # orthogonalize
  W = orth(W)
  evs = eigen(cov(X %*% W))
  evals = evs$values
  evecs = evs$vectors
  
  W = W %*% evecs
  Z = X %*% W

  if (showits)                                                     # Show last iteration
    cat(paste0(format(it), "...", "\n"))
  
  return(list(scores=Z, loadings=W, reconerr=reconerr, Xrecon=t(Xrecon)))
}

###############
### Example ###
###############

### Get data and run
# state.x77 is the data; various state demographics
X = scale(state.x77)
outEM = PCAEM(X=X, nComp=2, tol=1e-12, maxit=1000)
outEM

# Extract reconstructed values and loadings for comparison
Xrecon = outEM$Xrecon
loadingsEM = outEM$loadings
scoresEM = outEM$scores

# mean squared reconstruction error
mean((Xrecon-X)^2)  # outEM$reconerr/prod(dim(X))


### compare results to output from pcaMethods; note that signs for loadings/scores may be different
library(pcaMethods)
outpcam = pca(X, nPcs=2, method='svd', scale='none', center=F)
loadings_pcam = loadings(outpcam)
scores_pcam = scores(outpcam)

# compare loadings and scores
sum((abs(loadings_pcam)-abs(loadingsEM))^2)
abs(round(cbind(scores_pcam, scoresEM), 2))

# compare reconstructed data sets
Xrecon2 = scores_pcam %*% t(loadings_pcam)
mean((Xrecon2-X)^2)
mean(abs(Xrecon2-Xrecon))

# plots
plot(Xrecon2[,1], X[,1])
plot(Xrecon2[,2], X[,2])

plot(Xrecon[,1], Xrecon2[,1])
