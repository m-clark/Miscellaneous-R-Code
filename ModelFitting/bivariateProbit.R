### Stata's seems to be the primary audience concerned with these models, but I
### thought I'd play around with one here (I've never had reason to use a probit
### model in practice).  Stata examples come from the UCLA ATS website and the
### Stata manual so one can investigate the Stata result for comparison.



# standard probit --------------------------------------------------------

probitLL = function(beta, X, y){
  mu = X %*% beta
  
  # these produce identical results, but the second is the typical depiction
  # ll = sum(dbinom(
  #   y,
  #   size = 1,
  #   prob = pnorm(mu),
  #   log = T
  # ))     
  
  ll = sum(y * pnorm(mu, log = T) + (1 - y) * log(1 - pnorm(mu)))
  
  -ll
}

### Example 1

# Example at https://stats.idre.ucla.edu/stata/dae/probit-regression/

admit = haven::read_dta('https://stats.idre.ucla.edu/stat/stata/dae/binary.dta')

head(admit)

X = model.matrix(admit~ gre + gpa + factor(rank), admit)
y = admit$admit

init = rep(0, ncol(X))

optimResult = optim(
  fn  = probitLL,
  par = init,
  X = X,
  y = y,
  method = 'BFGS'
)

optimResult 


### Example 2: from Stata manual on standard probit 

# http://www.stata.com/manuals13/rprobit.pdf 
# "We have data on the make, weight, and mileage rating of 22 foreign and 52 
# domestic automobiles. We wish to fit a probit model explaining whether a car 
# is foreign based on its weight and mileage."

auto = haven::read_dta('http://www.stata-press.com/data/r13/auto.dta')

head(auto)

X = model.matrix(foreign~ weight + mpg, auto)
y = auto$foreign

init = rep(0, ncol(X))

optimResult = optim(
  fn  = probitLL,
  par = init,
  X = X,
  y = y
)

optimResult



# Bivariate probit --------------------------------------------------------

bivariateProbitLL = function(pars, X, y1, y2) {
  rho = pars[1]
  mu1 = X %*% pars[2:(ncol(X) + 1)]
  mu2 = X %*% pars[(ncol(X) + 2):length(pars)]
  q1 = ifelse(y1 == 1, 1,-1)
  q2 = ifelse(y2 == 1, 1,-1)
  
  require(mnormt)
  eta1 = q1 * mu1
  eta2 = q2 * mu2
  
  ll = matrix(NA, nrow = nrow(X))
  for (i in 1:length(ll)) {
    corr = q1[i] * q2[i] * rho
    corr = matrix(c(1, corr, corr, 1), 2)
    ll[i] = log(
      pmnorm(
        x = c(eta1[i], eta2[i]),
        mean = c(0, 0),
        varcov = corr
      )
    )
  }
  
  # the loop is probably clearer, and there is no difference in time, but here's a oneliner
  # ll = mapply(function(e1, e2, q1, q2) log(pmnorm(x=c(e1, e2), varcov = matrix(c(1,q1*q2*rho,q1*q2*rho,1),2))),
  #             eta1, eta2, q1, q2)
  
  -sum(ll)
}


### Example 3: from stata manual on bivariate probit
# "We wish to model the bivariate outcomes of whether children attend private 
# school and whether the head of the household voted for an increase in property
# tax based on the other covariates."

school = haven::read_dta('http://www.stata-press.com/data/r13/school.dta')

head(school)

X  = model.matrix(private ~ years + logptax + loginc, school)
y1 = school$private
y2 = school$vote

init = c(0, rep(0, ncol(X)*2))

# you'll probably get a warning or two, ignore; takes a couple seconds
optimResult = optim(
  fn  = bivariateProbitLL,
  par = init,
  X   = X,
  y1  = y1,
  y2  = y2,
  method = 'BFGS'
)


loglik = optimResult$value
rho = optimResult$par[1]
coefsPrivate = optimResult$par[2:(ncol(X) + 1)]
coefsVote = optimResult$par[(ncol(X) + 2):length(init)]
names(coefsPrivate) = names(coefsVote) = c('Int', 'years', 'logptax', 'loginc')

list(
  loglik  = loglik,
  rho     = rho,
  Private = coefsPrivate,
  Vote    = coefsVote
)

