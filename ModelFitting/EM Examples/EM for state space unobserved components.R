### The following regards chapter 11 in Statistical Modeling and Computation, 
### the first example for an unobserved components model.  The data regards 
### inflation based on the U.S. consumer price index (infl = 
### 400*log(cpi_t/cpi_{t-1})), from the second quarter of 1947 to the second 
### quarter of 2011.  You can acquire the data here 
### (http://www.maths.uq.edu.au/~kroese/statbook/Statespace/USCPI.csv) or in
### Datasets repo. just note that it has 2 mystery columns and one mystery row
### presumably supplied by Excel.  You can also get the CPI data yourself at
### http://www.bls.gov/cpi/ in a frustrating fashion, or in a much easier
### fashion here 
### http://research.stlouisfed.org/fred2/series/CPIAUCSL/downloaddata.


### For the following I use n instead of t or T because those are transpose and
### TRUE in R. The model is basically y = τ + ϵ, with ϵ ~ N(0, σ^2), 
### and τ = τ_{n-1} + υ_n with υ ~ N(0, ω^2).  Thus each y is
### associated with a latent variable that follows a random walk over time.
### ω^2 serves as a smoothing parameter, which itself may be estimated but
### which is fixed in the following. See the text for more details.

d = read.csv('../Datasets/USCPI.csv', header=F)
inflation = d[,1]
summary(inflation)


################################
### EM for state space model ###
################################
statespaceEM = function(params, y, omega2_0, omega2, tol=.00001, maxits=100, showits=T){     
  # Arguments are starting parameters (variance as 'sigma2'), data, tolerance,
  # maximum iterations, and whether to show iterations
  
  # Not really needed here, but would be a good idea generally to take advantage
  # of sparse representation for large data
  require(spam)
  
  # Starting points
  n = length(y)
  sigma2 = params$sigma2
  
  # Other initializations
  H = diag(n)
  for (i in 1:(ncol(H)-1)){
    H[i+1,i] = -1
  }
  
  Omega2 = as.spam(diag(omega2, n)); Omega2[1,1] = omega2_0
  H = as.spam(H)
  HinvOmega2H = t(H) %*% chol2inv(chol(Omega2)) %*% H      # tau ~ N(0, HinvOmmega2H^-1)
  
  
  it = 0
  converged = FALSE
  
  if (showits)                                             # Show iterations
    cat(paste("Iterations of EM:", "\n"))

  while ((!converged) & (it < maxits)) { 
    sigma2Old = sigma2[1]
    Sigma2invOld = diag(n)/sigma2Old

    K = HinvOmega2H + Sigma2invOld                         # E
    tau = solve(K, y/sigma2Old)                            # tau|y, sigma2_{n-1}, omega2 ~ N(0, K^-1)

    K_inv_tr = sum(1/eigen(K)$values)
    
    sigma2 = 1/n * (K_inv_tr + crossprod(y-tau))           # M
    
    converged = max(abs(sigma2 - sigma2Old)) <= tol
    
    # if showits true, & it =1 or divisible by 5 print message
    it = it + 1
    if (showits & it == 1 | it%%5 == 0)        
      cat(paste(format(it), "...", "\n", sep = ""))
  }
  
  Kfinal = HinvOmega2H + diag(n)/sigma2[1];
  taufinal = solve(K,(y/sigma2));
  out = list(sigma2=sigma2, tau=taufinal)
}

# debugonce(statespaceEM)
ssMod_1 = statespaceEM(params=data.frame(sigma2=var(inflation)), y=inflation, tol=1e-4, omega2_0=9, omega2=1^2)
ssMod_.5 = statespaceEM(params=data.frame(sigma2=var(inflation)), y=inflation, tol=1e-4, omega2_0=9, omega2=.5^2)
ssMod_1$sigma2
ssMod_.5$sigma2

library(lubridate)
series = ymd(paste0(rep(1947:2014, e=4),'-', c('01','04','07','10'), '-', '01'))
seriestext = series[1:length(inflation)]

# fig. 11.1 in the text
plot(seriestext, inflation, pch=19, cex=.5, col='gray50', ylim=c(-10,20), bty='n')
lines(seriestext, ssMod_1$tau, col="#FF5500")
lines(seriestext, ssMod_.5$tau, col="skyblue3")
