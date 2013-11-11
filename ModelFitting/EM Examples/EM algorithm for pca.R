PCAEM = function(X, nfactors=2, classic=F, tol=.00001, maxits=100, showits=T){       # tol = tolerance
  require(pracma) # for orthonormal basis of W
  
  #starting points
  N = nrow(X)
  D = ncol(X)
  L = nfactors
  Xt = t(X)
  Z = t(replicate(L, rnorm(N)))
  W = replicate(L, rnorm(D))
  
  it <- 0
  converged <- FALSE
  
  
  if (showits)                                                            # Show iterations
    cat(paste("Iterations of EM:", "\n"))
  while ((!converged) & (it < maxits)) {                                  # while no convergence and we haven't reached our max iterations do this stuff
    Z.old <- Z                                                            # create 'old' values for comparison
    Z = solve(t(W)%*%W) %*% crossprod(W, Xt)                                     # E
    W = Xt%*%t(Z) %*% solve(tcrossprod(Z))                                # M

    it <- it + 1
    if (showits & (it == 1 | it%%5 == 0))
      cat(paste(format(it), "...", "\n", sep = ""))
    converged <- max(abs(Z.old-Z)) <= tol
  }
  # calculate reconstruction error
  Xrecon = W %*% Z
  reconerr = sum((Xrecon-t(X))^2)
  
  # orthogonalize
  W <- orth(W)
  evs <- eigen(cov(X %*% W))
  evals = evs$values
  evecs = evs$vectors
  
  W <- W %*% evecs
  Z <- X %*% W

  if (showits)                                                            # Show last iteration
    cat(paste0(format(it), "...", "\n"))
  
  return(list(scores=Z, loadings=W,  reconerr = reconerr, Xrecon=t(Xrecon)))
}

# faithful  state.x77
X = scale(state.x77)
outEM = PCAEM(X=X, nfactors=4, tol=1e-12, maxit=1000)
outEM
cor(outEM$scores) 
Xrecon = outEM$Xrecon
mean((Xrecon-X)^2)  

plot(Xrecon[,1], X[,1])
plot(Xrecon[,2], X[,2])


library(pcaMethods)
outpcam = pca(X, nPcs=4, method='svd')
list(loadings(outpcam), outEM$loadings)
cbind(scores(outpcam), outEM$scores)
Xrecon2 = loadings(outpcam) %*% t(scores(outpcam))
Xrecon2 = t(Xrecon2)
mean((Xrecon2-X)^2)
mean(abs(Xrecon2-Xrecon))

plot(Xrecon2[,1], X[,1])
plot(Xrecon2[,2], X[,2])

plot(Xrecon[,1], Xrecon2[,1])