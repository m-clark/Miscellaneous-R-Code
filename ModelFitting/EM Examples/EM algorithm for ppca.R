PPCAEM = function(X, nfactors=2, classic=F, tol=.00001, maxits=100, showits=T){       # tol = tolerance
  require(pracma) # for orthonormal basis of W
  
  #starting points
  N = nrow(X)
  D = ncol(X)
  L = nfactors
  Xt = t(X)
  Z = t(replicate(L, rnorm(N)))
  
  S = (1/N) * Xt%*%X
  evals = eigen(S)$values
  evecs = eigen(S)$vectors

  V = evecs[,1:L]
  Lambda = diag(evals[1:L])
  sigma2 = (1/(D-L)) * sum(evals[(L+1):D])
  
  W = V %*% chol(Lambda-sigma2*diag(L)) %*% diag(L)
#   Psi = sigma2*diag(D)  # D x D 
#   C = W%*%t(W) + Psi
  
  
  it <- 0
  converged <- FALSE
  ll = 0
  
  if (showits)                                                            # Show iterations
    cat(paste("Iterations of EM:", "\n"))
  while ((!converged) & (it < maxits)) {                                  # while no convergence and we haven't reached our max iterations do this stuff
    Z.old <- Z                                                            # create 'old' values for comparison
    ll.old = ll
    Z = solve(t(W)%*%W) %*% t(W)%*%Xt
    Sigmahat = cov(Z)
#     W = V %*% chol(Lambda-sigma2*diag(L)) %*% diag(L)

    W = Xt%*%t(Z) %*% solve(tcrossprod(Z))
    W = orth(W)
    Psi = sigma2*diag(D)  # D x D 
    C = W%*%t(W) + Psi
    ll = (-N/2) * log(det(C)) + sum(diag(solve(C)%*%S))
    

    Xrecon = W %*% Z
    
    it <- it + 1
    if (showits & (it == 1 | it%%5 == 0))
      cat(paste(format(it), "...", "\n", sep = ""))
    converged <- max(abs(Z.old-Z)) <= tol
  }
  reconerr = sum((Xrecon-t(X))^2)
  
  if (showits)                                                            # Show last iteration
    cat(paste0(format(it), "...", "\n"))
  
  return(list(scores=t(Z), loadings=W,  reconerr = reconerr))
}


X = scale(state.x77)
outEM = PPCAEM(X=X, nfactors=2, tol=1e-8, maxit=100)
outEM  
Xrecon = outEM$loadings %*% t(outEM$scores )
Xrecon = t(Xrecon)
mean((Xrecon-X)^2)

plot(Xrecon2[,1], X[,1])
plot(Xrecon2[,2], X[,2])


library(pcaMethods)
outpcam = ppca(X, nPcs=2)
list(loadings(outpcam), outEM$loadings)
cbind(scores(outpcam), outEM$scores)
Xrecon2 = loadings(outpcam) %*% t(scores(outpcam))
Xrecon2 = t(Xrecon2)
mean((Xrecon2-X)^2)
mean(abs(Xrecon2-Xrecon))