### topic models
### IN PROGRESS AND INCOMPLETE

library(topicmodels)
data('AssociatedPress')

# The example code in the manual/online proposes a 'bag of words' representation
# instead of working with a document term matrix as would perhaps be the usual
# case.  I first attempt to see if the code as posited works.
ap = as.matrix(AssociatedPress)[1:100,]  # first 100 docs


# convert the DTM to lists of word IDs as a programming exercise
library(plyr)
apBoW = alply(ap, 1, function(wordcounts) rep(colnames(ap), wordcounts))
apBoWIDs = lapply(apBoW, function(wordlist) which(colnames(ap) %in% wordlist)[factor(wordlist)])

# examine, because yes that was only two lines of code with no loop
rbind(apBoW[[12]], apBoWIDs[[12]])[,1:20]

# alternatively one could extract directly from the simple triplet matrix DTM
# apSTM = data.frame(doc=AssociatedPress$i, word=AssociatedPress$j, count=AssociatedPress$v)
# apSTM = dplyr::filter(apSTM, doc %in% 1:100)
# w = rep(apSTM$word, apSTM$count)

w = unlist(apBoWIDs)
V = ncol(ap)
M = nrow(ap)
N = length(w)
doc = rep(1:M, sapply(apBoWIDs, length))

K = 10
alpha = rep(1/K, K)
beta = colSums(ap)/sum(ap)
beta = rep(1/V, V)

standat = list(K=K, V=V, M=M, N=N, w=w, doc=doc, alpha=alpha, beta=beta)

### Stan manual for LDA
stanLDA <- '
data {
  int<lower=2> K;                                // num topics
  int<lower=2> V;                                // num words
  int<lower=1> M;                                // num docs
  int<lower=1> N;                                // total word instances
  int<lower=1,upper=V> w[N];                     // word n
  int<lower=1,upper=M> doc[N];                   // doc ID for word n
  vector<lower=0>[K] alpha;                      // topic prior
  vector<lower=0>[V] beta;                       // word prior
}

parameters {
  simplex[K] theta[M];                           // topic dist for doc m
  simplex[V] phi[K];                             // word dist for topic k
}

model {
  for (m in 1:M)
    theta[m] ~ dirichlet(alpha);                 // prior
  for (k in 1:K)
    phi[k] ~ dirichlet(beta);                    // prior
  for (n in 1:N) {
    real gamma[K];
    for (k in 1:K)
      gamma[k] <- log(theta[doc[n],k]) + log(phi[k,w[n]]);
    increment_log_prob(log_sum_exp(gamma));      // likelihood
  }
}
'

### test run
library(rstan)
library(parallel)
cl = makeCluster(2)
clusterEvalQ(cl, library(rstan))
clusterExport(cl, c('standat', 'stanLDA'))

p = proc.time()
test = parSapply(cl, 1:2, function(i) stan(model_code=stanLDA, data=standat, chains=1, chain_id=i))
proc.time() - p

stopCluster(cl)

testOut = sflist2stanfit(test)

print(testOut, par='theta')

thetas = extract(testOut)$theta
phi = extract(testOut)$phi










### Stan manual for CTM
stanCTM <- '
data {
  int<lower=2> K;                                // num topics
  int<lower=2> V;                                // num words
  int<lower=1> M;                                // num docs
  int<lower=1> N;                                // total word instances
  int<lower=1,upper=V> w[N];                     // word n
  int<lower=1,upper=M> doc[N];                   // doc ID for word n
  vector<lower=0>[V] beta;                       // word prior
}

parameters {
  vector[K] mu;                                  // topic mean
  corr_matrix[K] Omega;                          // correlation matrix
  vector<lower=0>[K] sigma;                      // scales
  vector[K] eta[M];                              // logit topic dist for doc m
  simplex[V] phi[K];                             // word dist for topic k
}

transformed parameters {
  simplex[K] theta[M];                           // simplex topic dist for doc m
  cov_matrix[K] Sigma;                           // covariance matrix

  for (m in 1:M)
    theta[m] <- softmax(eta[m]);
  
  for (m in 1:K) {
    Sigma[m,m] <- sigma[m] * sigma[m] * Omega[m,m];
    for (n in (m+1):K) {
      Sigma[m,n] <- sigma[m] * sigma[n] * Omega[m,n];
      Sigma[n,m] <- Sigma[m,n];
    }
  }
}

model {
  // priors
  for (k in 1:K)
    phi[k] ~ dirichlet(beta);
  mu ~ normal(0,5);
  Omega ~ lkj_corr(2.0);
  sigma ~ cauchy(0,5);
  // topic distribution for docs
  for (m in 1:M)
    eta[m] ~ multi_normal(mu,Sigma);
  // token probabilities
  for (n in 1:N) {
  real gamma[K];
    for (k in 1:K)
      gamma[k] <- log(theta[doc[n],k]) + log(phi[k,w[n]]);
      increment_log_prob(log_sum_exp(gamma));    // likelihood
    }
  }
'
