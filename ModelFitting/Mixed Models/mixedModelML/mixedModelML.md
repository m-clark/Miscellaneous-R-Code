# Mixed Model Estimation via Maximum Likelihood




## Introduction
This is an example based on Wood, 2006, chapter 6 in particular.  It assumes familiarity with standard regression from a matrix perspective. 

## Model Formulation

We can start with a standard linear model expressed as follows:

$$\mathbf{y} = \mathbf{Xb} + \mathbf{\epsilon} $$

Here $\mathbf{y}$ is the target variable, $\mathbf{X}$ is a model matrix (first column representing the intercept, the rest are the covariates of interest), $\mathbf{b}$ are the coefficients, and error $\mathbf{\epsilon}$. Note that beyond this point I won't use bold to indicate vectors/matrices, nor subscripts for every *i*<sup>th</sup> observation. Let's just assume we are in a normal data situation involving more than one observation, a univariate vector target variable (y), a matrix of predictor variables (X) etc.

For a mixed model with a single random effect for some grouping factor (e.g. students within schools), this extends to:

$$y = Xb + Zg + \epsilon$$

Where Z is an indicator matrix pertaining to the grouping structure (sometimes referred to as dummy coding or one-hot encoding). Consider a factor z representing group/cluster membership, this would convert z to the following:


------------------
 z   ZA   ZB   ZC 
--- ---- ---- ----
 A   1    0    0  

 A   1    0    0  

 B   0    1    0  

 B   0    1    0  

 C   0    0    1  

 C   0    0    1  
------------------

The coefficients $g$ are the random effects, assumed $\mathcal{N}(0,\tau)$, and while we are often interested in them, they do not have to be estimated directly.

$$y = Xb + Zg + \epsilon \\
g \sim \mathcal{N}(0, \psi_\theta) \\
\epsilon \sim \mathcal{N}(0, \Lambda\sigma^2)$$

In this depiction $\psi_\theta$ can reflect some more interesting dependencies, but in the simple case of a random intercepts model it can be a single variance estimate $\tau^2$. $\Lambda$ can be used to model residual covariance but often is just the identity matrix, with the underlying assumption of constant variance $\sigma^2$ across observations.  

We can combine the random and residuals into a single construct reflecting the covariance structure of the observations:

$$ e = Zg + \epsilon $$

This makes $\mathbf{e}$ a multivariate vector with mean 0 and covariance:

$$Z\psi_{\theta}Z^\intercal + I\sigma^2$$

This puts us back to a standard linear model:

$$ y = Xb + e, \\
e \sim \mathcal{N}(0, \Sigma_\theta\sigma^2)$$


## Maximum Likelihood Estimation

Given where we are now, we can proceed to estimate the mixed model. For this we'll use the sleepstudy data from lme4. The data has reaction times for 18 individuals over 10 days each (see the help file for the sleepstudy object for more details).

### Data

```r
data(sleepstudy, package='lme4')
X = model.matrix(~Days, sleepstudy)
Z = model.matrix(~factor(sleepstudy$Subject)-1)
y = sleepstudy$Reaction
```

### ML function
The following is based on the code in Wood (6.2.2), with a couple modifications for consistent nomenclature. $\theta$ represents the vector of parameters we wish to estimate. The (square root of the) variances will be estimated on the log scale. In Wood, he simply extracts the 'fixed effects' for the intercept and days effects using lm (6.2.3).


```r
llMixed = function(y, X, Z, theta){
  tau = exp(theta[1])
  sigma = exp(theta[2])
  n = length(y)
  
  # evaluate cov mat for y
  e = tcrossprod(Z)*tau^2 + diag(n)*sigma^2
  L = chol(e)  # L'L = e
  
  # transform dependent linear model to independent
  y = backsolve(L, y, transpose=TRUE)
  X = backsolve(L, X, transpose=TRUE)
  b = coef(lm(y~X-1))
  LP = X %*% b
  
  ll = -n/2*log(2*pi) -sum(log(diag(L))) - crossprod(y-LP)/2
  -ll
}
```


Here is an alternative function using a multivariate approach that doesn't use the transformation to independent, and might provide additional perspective..


```r
llMixedMV = function(y, X, Z, theta){
  tau = exp(theta[1])
  sigma = exp(theta[2])
  n = length(y)
  
  # evaluate cov mat for y
  e = tcrossprod(Z)*tau^2 + diag(n)*sigma^2

  b = coef(lm.fit(X, y))
  mu = X %*% b

  ll = -mvtnorm::dmvnorm(y, mu, e, log=T)
}
```



### Results

We'll use the optim function for estimation.  A slight change to tolerance is included to get a closer estimate to lme4.


```r
paramInit = c(0, 0)
names(paramInit) = c('tau', 'sigma')

modelResults = optim(llMixed, X=X, y=y, Z=Z, par=paramInit, control=list(reltol=1e-10))
modelResultsMV = optim(llMixedMV, X=X, y=y, Z=Z, par=paramInit, control=list(reltol=1e-10))

rbind(c(exp(modelResults$par), logLik = modelResults$value, coef(lm(y~X-1))),
      c(exp(modelResultsMV$par), logLik = modelResultsMV$value, coef(lm(y~X-1)))) %>% 
  round(2)
```

```
##        tau sigma logLik X(Intercept) XDays
## [1,] 36.02  30.9 897.04       251.41 10.47
## [2,] 36.02  30.9 897.04       251.41 10.47
```

As we can see, both formulations produce identical results. We can now compare those results to the lme4 output for the same model.


```r
library(lme4)
lmeMod = lmer(Reaction ~ Days + (1|Subject), sleepstudy, REML=FALSE)
lmeMod
```

```
## Linear mixed model fit by maximum likelihood  ['lmerMod']
## Formula: Reaction ~ Days + (1 | Subject)
##    Data: sleepstudy
##       AIC       BIC    logLik  deviance  df.resid 
## 1802.0786 1814.8505 -897.0393 1794.0786       176 
## Random effects:
##  Groups   Name        Std.Dev.
##  Subject  (Intercept) 36.01   
##  Residual             30.90   
## Number of obs: 180, groups:  Subject, 18
## Fixed Effects:
## (Intercept)         Days  
##      251.41        10.47
```

We can predict the random effects (Wood, 6.2.4), and after doing so compare to the lme4 estimates.


```r
tau = exp(modelResults$par)[1]
tausq = tau^2
sigma = exp(modelResults$par)[2]
sigmasq = sigma^2
Sigma = tcrossprod(Z)*tausq/sigmasq + diag(length(y))
ranefEstimated = tausq*t(Z)%*%solve(Sigma) %*% resid(lm(y~X-1))/sigmasq
data.frame(ranefEstimated, lme4 = ranef(lmeMod)$Subject[[1]]) %>% round(2)
```

```
##                               ranefEstimated   lme4
## factor(sleepstudy$Subject)308          40.64  40.64
## factor(sleepstudy$Subject)309         -77.57 -77.57
## factor(sleepstudy$Subject)310         -62.88 -62.88
## factor(sleepstudy$Subject)330           4.39   4.39
## factor(sleepstudy$Subject)331          10.18  10.18
## factor(sleepstudy$Subject)332           8.19   8.19
## factor(sleepstudy$Subject)333          16.44  16.44
## factor(sleepstudy$Subject)334          -2.99  -2.99
## factor(sleepstudy$Subject)335         -45.12 -45.12
## factor(sleepstudy$Subject)337          71.92  71.92
## factor(sleepstudy$Subject)349         -21.12 -21.12
## factor(sleepstudy$Subject)350          14.06  14.06
## factor(sleepstudy$Subject)351          -7.83  -7.83
## factor(sleepstudy$Subject)352          36.25  36.25
## factor(sleepstudy$Subject)369           7.01   7.01
## factor(sleepstudy$Subject)370          -6.34  -6.34
## factor(sleepstudy$Subject)371          -3.28  -3.28
## factor(sleepstudy$Subject)372          18.05  18.05
```


I'd like to at some point demonstrate some concepts from section 6.6 in Wood, but that will have to wait for now.
