# markov chains and hidden markov models




# markovchain package demo ------------------------------------------------

library(markovchain)
A = matrix(c(.7,.3,.9,.1), nrow=2, byrow=T)
dtmcA = new('markovchain', transitionMatrix=A,
             states=c('a','b'),
             name='MarkovChain A')
dtmcA
plot(dtmcA)
transitionProbability(dtmcA, 'b', 'b')
initialState = c(0,1)
steps = 4
finalState = initialState*dtmcA^steps #using power operator
finalState


steadyStates(dtmcA)


observed_states = sample(c('a', 'b'), 50, c(.7, .3), replace=T)
createSequenceMatrix(observed_states)
markovchainFit(observed_states)



# Create data -------------------------------------------------------------

# a recursive function to take a matrix power
mat_power<- function(M, N){
  if (N==1) return(M)
  
  M %*% mat_power(M, N-1)
}

# example
test.mat = matrix(rep(2,4), nrow=2)
mat_power(test.mat, 2)

# transition matrix
A = matrix(c(.7,.3,.4,.6), nrow=2, byrow=T)
mat_power(A, 10)


# a function to create a sequence
createSequence = function(states, len, tmat) {
  # states: number of states
  # len: length of sequence
  # tmat: the transition matrix
  states_numeric = length(unique(states))
  out = numeric(len)
  out[1] = sample(states_numeric, 1, prob=colMeans(tmat)) # initial state
  
  for (i in 2:len){
    out[i] = sample(states_numeric, 1, prob=tmat[out[i-1],])
  }
  states[out]
}



# Two state demo ----------------------------------------------------------

# Note that a notably long sequence is needed to get close to recovering the
# true transition matrix
A = matrix(c(.7,.3,.9,.1), nrow=2, byrow=T)
observed_states = createSequence(c('a', 'b'), 5000, tmat=A)
createSequenceMatrix(observed_states)
prop.table(createSequenceMatrix(observed_states), 1)

markovchainFit(observed_states)
res = markovchainFit(observed_states)

# log likelihood
sum(createSequenceMatrix(observed_states) * log(res$estimate@transitionMatrix))



# Three state demo --------------------------------------------------------

A = matrix(c(.7,.2, .1,
             .2, .4, .4,
             .05,.05, .9), nrow=3, byrow=T)

observed_states = createSequence(c('a', 'b', 'c'), 500, tmat=A)
createSequenceMatrix(observed_states)
prop.table(createSequenceMatrix(observed_states), 1)
markovchainFit(observed_states)



# Fit a Markov Model ------------------------------------------------------

# Now we create a function to calculate the (negative) log likeihood

markov_model <- function(par, x) {
  # par should be the c(A) of tran probabilities A
  nstates = length(unique(x))
  
  # create transition matrix
  par = matrix(par, ncol=nstates)
  par = t(apply(par, 1, function(x) x/sum(x)))
  
  # create seq matrix
  seqMat = table(x[-length(x)], x[-1])
  
  # calculate log likelihood
  ll = sum(seqMat*log(par))
  -ll
}

A = matrix(c(.7,.2, .1,
             .40, .2, .40,
             .1,.15,.75), nrow=3, byrow=T)
observed_states = createSequence(c('a', 'b', 'c'), 1000, tmat=A)


# note that initial state values will be transformed to rowsum to one, so the 
# specific initial values don't matter (i.e. don't have to be probabilities). 
# With the basic optim approach, sometimes log(0) will occur and produce
# warning. Can be ignored, or use LFBGS as below.

initpar = rep(1, 9)
test = optim(initpar, markov_model, x=observed_states, method='BFGS',
             control=list(reltol=1e-12))

# get estimates on prob scale
estmat = matrix(test$par, ncol=3)
estmat = t(apply(estmat, 1, function(x) x/sum(x)))

# compare with markov chain package
compare_result = markovchainFit(observed_states)

# compare log likelihood
c(-test$value, compare_result$logLikelihood)

# compare estimated transition matrix
list(`Estimated via optim`= estmat, 
     `markovchain Package`= compare_result$estimate@transitionMatrix,
     `Analytical Solution`= prop.table(table(observed_states[-length(observed_states)], 
                                             observed_states[-1]), 1)) %>% 
  lapply(round, 3)

# plot
plot(new('markovchain', transitionMatrix=estmat,
    states=c('a','b', 'c'),
    name='Estiamted Markov Chain'))


# if you don't want warnings due to zeros; see also constrOptim
# test = optim(initpar, markov_model, x=observed_states, method='L-BFGS', 
#              lower=rep(1e-20, length(initpar)), 
#              control=list(pgtol=1e-12))
