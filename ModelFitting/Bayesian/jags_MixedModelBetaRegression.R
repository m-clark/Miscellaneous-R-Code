### The following is an attempt to reproduce the Stan approach to a mixed model 
### with a beta distribution for the response response in
### rstan_MixedModelBetaRegression.R.  For more detail about what the JAGS code
### is doing, see the Stan code in that file.  See the helpfile for
### GasolineYield for information on the data. 

### A couple of changes to note. JAGS doesn't have the Cauchy distribution, so I
### show the approach in the code outlined in Gelman 2006.  Instead of that, one
### can also use an approach with the t-distribution similarly truncated. Unlike
### the Stan code, the variance components and phi priors are all half.cauchy(0,
### 5).  I'll have to do a little more playing with it but it seemed initially
### that JAGS took a lot more iterations to start looking as healthy as the Stan
### results in terms of effective sample size (though this is the case for many
### models), having particular trouble with the random intercepts and variance
### estimate (though Stan had issues there too).  But just for comparison's
### sake, it will reproduce Stan parameter results using the same amount of
### iterations, so I leave it shorter.  Keep in mind also we're not dealing with
### a lot of data here.

##################
### Data setup ###
##################
data('GasolineYield', package='betareg')
jagsdat = list('yield'=GasolineYield$yield, 'tempCen'=c(scale(GasolineYield$temp, scale=F)), 
               'N'=nrow(GasolineYield), 'id'=as.numeric(GasolineYield$batch), L=10)


#######################
### Jags Model Code ###
#######################

# write to file
sink('betaMixedModeljags.txt')
cat(
  'model {
  for (n in 1:N){
    logit(mu[n]) <- (Intercept + IntRE[id[n]]) + (betaTemp + SlopeRE[id[n]])*tempCen[n]
    A[n] <- mu[n]*phi
    B[n] <- (1.0-mu[n])*phi
    yield[n] ~ dbeta(A[n], B[n])
  }

  Intercept ~ dnorm(0, 1/10^2)
  betaTemp ~ dnorm(0, 1)
  
  for (l in 1:L){
    IntRE[l] <- gammaIntercept[l]*sd_int
    SlopeRE[l] <- gammaTemp[l]*sd_beta

    gammaIntercept[l] ~ dnorm(0, 1)
    gammaTemp[l] ~ dnorm(0, 1)
  }

  # Half-cauchy as in Gelman 2006
  # If scale parameter in cauchy is 5, precision of z = 1/5^2 = 0.04
  # sigma int
  #   sd_int <- zInt/sqrt(chSqInt)                                                  # prior for sigma; cauchy = normal/sqrt(chi^2)
  #   zInt ~ dnorm(0, .04)I(0,)
  #   chSqInt ~ dgamma(0.5, 0.5)                                                    # chi^2 with 1 d.f.
  sd_int ~ dt(0, .04, 1)I(0,)

  # sigma beta
  #   sd_beta <- zbeta/sqrt(chSqbeta)                                               
  #   zbeta ~ dnorm(0, .04)I(0,)
  #   chSqbeta ~ dgamma(0.5, 0.5)                                                   
  sd_beta ~ dt(0, .04, 1)I(0,)

  # phi
  #   phi <- zphi/sqrt(chSqphi)                                                     
  #   zphi ~ dnorm(0, .04)I(0,)
  #   chSqphi ~ dgamma(0.5, 0.5)       
  phi ~ dt(0, .04, 1)I(0,)
}'
)
sink()



#####################
### Run the model ###
#####################
library(rjags)
mixedbeta0 <- jags.model(file='betaMixedModeljags.txt', data=jagsdat, # inits=inits
                           n.chains=4, n.adapt=2000)

# update(mixedbetamod, 10000)

mixedbeta = coda.samples(mixedbeta0, c('Intercept','betaTemp','IntRE', 'SlopeRE', 
                                       'phi', 'sd_int', 'sd_beta'), n.iter=200000, 
                         thin=200, n.chains=4)
summary(mixedbeta)

### cleaner view of coefficients of interest
round(summary(mixedbeta)[[1]], 4)

### effective sample size
effectiveSize(mixedbeta)

### trace and density plots
plot(mixedbeta)
