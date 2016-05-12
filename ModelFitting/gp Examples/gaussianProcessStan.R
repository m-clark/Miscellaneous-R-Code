# Data and parameter setup ------------------------------------------------

# Data
set.seed(1234)
N = 20
Ntest = 200
x = rnorm(N, sd=1)
y = scale(sin(x) + rnorm(N, sd=.1))[,1]
xtest = seq(min(x)-1, max(x)+1, l=Ntest)
plot(x,y, pch=19, col='#ff5500')

# parameters
eta_sq = 1
rho_sq = 1
sigma_sq = .1

# Covariance function same as implemented in the Stan code.
Kfn <- function (x, eta_sq, rho_sq, sigma_sq) {
  N = length(x)
  Sigma = matrix(NA, N, N)
  
  # off diag elements
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      Sigma[i,j] <- eta_sq * exp(-rho_sq * (x[i] - x[j])^2);
      Sigma[j,i] <- Sigma[i,j];
    }
  }
  
  # diagonal elements
  for (k in 1:N)
    Sigma[k,k] <- eta_sq + sigma_sq; # + jitter
  Sigma
}



# Vis prior functions -----------------------------------------------------
xinit = seq(-5,5,.2)
xprior = MASS::mvrnorm(3, 
                       mu=rep(0, length(xinit)), 
                       Sigma=Kfn(x=xinit,
                                 eta_sq = eta_sq, 
                                 rho_sq = rho_sq, 
                                 sigma_sq = sigma_sq))


library(reshape2)
gdat = melt(data.frame(x=xinit, y=t(xprior)), id='x')

library(ggvis)
gdat %>% 
  ggvis(~x, ~value) %>% 
  group_by(variable) %>% 
  layer_paths(strokeOpacity:=.5) %>% 
  add_axis('x', grid=F) %>% 
  add_axis('y', grid=F)



# Stan model code ---------------------------------------------------------
# models/covariance functions available
# gpStanModelCode_generalizedSquaredExponential.stan
# gpStanModelCode_gammaExponential.stan
# gpStanModelCode_rationalQuadratic.stan

gp = 'ModelFitting/gp Examples/gpStanModelCode_generalizedSquaredExponential.stan'

# Compile Check -----------------------------------------------------------

standata = list(N=N, x=x, y=y, xtest=xtest, Ntest=200)


library(rstan)
fit0 = stan(data=standata, file = gp, iter = 1, chains=1)



# Main Run ----------------------------------------------------------------

iterations = 12000
wu = 2000
th = 20
chains = 4


# With N = 20 Ntest = 200 takes about 2 min for gen squared expo, 10+min for 300

p = proc.time()
fit = stan(data=standata, file=gp, iter = iterations, warmup = wu, thin=th, 
           chains=chains, fit = fit0, cores=chains)
(proc.time() - p)/3600



# Summarize and Vis -------------------------------------------------------

# takes a bit to print
# print(fit, par=c('eta_sq','rho_sq','sigma_sq'))

# library(shinyStan)
# launch_shinystan(fit)


# Extract and visualize posterior predictive draws
yRep = extract(fit, 'yRep')$yRep

gdat = data.frame(x,y)
gdat2 = melt(data.frame(x = sort(xtest), y=t(yRep[sample(2000, 3),])), id='x')

gdat2 %>% 
  ggvis(~x, ~value) %>% 
  group_by(variable) %>% 
  layer_paths(strokeOpacity:=.25) %>% 
  layer_points(x=~x, y=~y, fill:='#ff5500', data=gdat) %>% 
  add_axis('x', grid=F) %>% 
  add_axis('y', grid=F)  


# Visualize fit
yRepMean = get_posterior_mean(fit, 'yRep')[,5]
quantiles = data.frame(t(apply(yRep, 2, quantile, p=c(.025,.975)))); colnames(quantiles) = c('ll','ul')

gdat3 = data.frame(x = sort(xtest), y=yRepMean, quantiles)
ttle = stringr::str_extract(gp, "(?<=_)(.*)(?=\\.)")

library(kernlab)
comparisonModel = gausspr(x, y, kernel=rbfdot(sigma=1), var=.1, scaled=F, tol=1e-5)
comparisonY = predict(comparisonModel, newdata=xtest)

library(mgcv)
comparisonModel2 = gam(y~s(x, bs='gp', m=c(2,2☻,2)))
comparisonY2 = predict(comparisonModel2, newdata=data.frame(x=xtest))

gdat3 %>% 
  ggvis(~x, ~y) %>% 
  layer_ribbons(y=~ll, y2=~ul, fillOpacity:=.1) %>% 
  layer_paths(strokeOpacity:=.5, stroke:='blue') %>%
  layer_paths(y=~comparisonY, strokeOpacity:=.5, stroke:='red') %>%
  layer_paths(y=~comparisonY2, strokeOpacity:=.5, stroke:='green') %>%
  layer_points(x=~x, y=~y, fill:='#ff5500', fillOpacity:=.5, data=gdat) %>% 
  add_axis('x', grid=F) %>% 
  add_axis("x", orient = "top", ticks = 0, title = ttle,
           properties = axis_props(
             axis = list(stroke = "white"),
             labels = list(fontSize = 0))) %>% 
  add_axis('y', grid=F)  


