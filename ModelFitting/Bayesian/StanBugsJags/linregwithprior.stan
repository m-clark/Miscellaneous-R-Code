data {                      // Data block; declarations only
  int<lower=0> N;           // Sample size                         
  int<lower=0> k;           // Dimension of model matrix
  matrix [N, k] X;          // Model Matrix
  vector[N] y;              // response
}

/* transformed data {       // Transformed data block; declarations and statements. None needed here.
 }
*/

parameters {                // Parameters block; declarations only
  vector[k] beta;           // coefficient vector
  real<lower=0> sigma;      // error scale
}

transformed parameters {    // Transformed parameters block; declarations and statements.
}

model {                     // Model block; declarations and statements.
  vector[N] mu;
  mu <- X * beta;           // creation of linear predictor

  // priors
  beta ~ normal(0, 10);
  sigma ~ cauchy(0, 5);     // With sigma bounded at 0, this is half-cauchy 

  // likelihood
  y ~ normal(mu, sigma);
}

generated quantities {      // Generated quantities block; declarations and statements.
  real rss;                
  real totalss;
  real R2;                  // Calculate Rsq as a demonstration
  vector[N] mu;
  
  mu <- X * beta;
  rss <- dot_self(mu-y);
  totalss <- dot_self(y-mean(y));
  R2 <- 1 - rss/totalss;
}