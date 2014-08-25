# From Translating Probability Density Functions: From R to BUGS and Back Again
# by David S. LeBauer, Michael C. Dietze, Benjamin M. Bolker

r2bugs.distributions <- function(priors, direction = 'r2bugs') {
  priors$distn <- as.character(priors$distn)
  priors$parama <- as.numeric(priors$parama)
  priors$paramb <- as.numeric(priors$paramb)
  ## index dataframe according to distribution
  norm <- priors$distn %in% c('norm', 'lnorm') # these have same transform
  weib <- grepl("weib", priors$distn) # matches r and bugs version
  gamma <- priors$distn == 'gamma'
  chsq <- grepl("chisq", priors$distn) # matches r and bugs version
  bin <- priors$distn %in% c('binom', 'bin') # matches r and bugs version
  nbin <- priors$distn %in% c('nbinom', 'negbin') # matches r and bugs version
  ## Normal, log-Normal: Convert sd to precision
  exponent <- ifelse(direction == "r2bugs", -2, -0.5)
  priors$paramb[norm] <- priors$paramb[norm] ^ exponent
  ## Weibull
  if(direction == 'r2bugs'){
    ## Convert R parameter b to BUGS parameter lambda by l = (1/b)^a
    priors$paramb[weib] <- (1 / priors$paramb[weib]) ^ priors$parama[weib]
  } else if (direction == 'bugs2r') {
    ## Convert BUGS parameter lambda to BUGS parameter b by b = l^(-1/a)
    priors$paramb[weib] <- priors$paramb[weib] ^ (- 1 / priors$parama[weib] )
  }
  ## Reverse parameter order for binomial and negative binomial
  priors[bin | nbin, c('parama', 'paramb')] <-
    priors[bin | nbin, c('paramb', 'parama')]
  ## Translate distribution names
  if(direction == "r2bugs"){
    priors$distn[weib] <- "weib"
    priors$distn[chsq] <- "chisqr"
    priors$distn[bin] <- "bin"
    priors$distn[nbin] <- "negbin"
  } else if(direction == "bugs2r"){
    priors$distn[weib] <- "weibull"
    priors$distn[chsq] <- "chisq"
    priors$distn[bin] <- "binom"
    priors$distn[nbin] <- "nbinom"
  }
  return(priors)
}

### Example ###
r.distn <- data.frame(distn = "norm", parama = 10, paramb = 2)
r2bugs.distributions(r.distn)