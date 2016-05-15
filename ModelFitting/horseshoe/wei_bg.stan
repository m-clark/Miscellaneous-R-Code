/*  Variable naming:
 obs       = observed
 cen       = (right) censored
 N         = number of samples
 M         = number of covariates
 bg        = established risk (or protective) factors
 tau       = scale parameter
*/
// Tomi Peltola, tomi.peltola@aalto.fi

functions {
  vector sqrt_vec(vector x) {
    vector[dims(x)[1]] res;

    for (m in 1:dims(x)[1]){
      res[m] <- sqrt(x[m]);
    }

    return res;
  }

  vector bg_prior_lp(real r_global, vector r_local) {
    r_global ~ normal(0.0, 10.0);
    r_local ~ inv_chi_square(1.0);

    return r_global * sqrt_vec(r_local);
  }
}

data {
  int<lower=0> Nobs;
  int<lower=0> Ncen;
  int<lower=0> M_bg;
  vector[Nobs] yobs;
  vector[Ncen] ycen;
  matrix[Nobs, M_bg] Xobs_bg;
  matrix[Ncen, M_bg] Xcen_bg;
}

transformed data {
  real<lower=0> tau_mu;
  real<lower=0> tau_al;

  tau_mu <- 10.0;
  tau_al <- 10.0;
}

parameters {
  real<lower=0> tau_s_bg_raw;
  vector<lower=0>[M_bg] tau_bg_raw;

  real alpha_raw;
  vector[M_bg] beta_bg_raw;

  real mu;
}

transformed parameters {
  vector[M_bg] beta_bg;
  real alpha;

  beta_bg <- bg_prior_lp(tau_s_bg_raw, tau_bg_raw) .* beta_bg_raw;
  alpha <- exp(tau_al * alpha_raw);
}

model {
  yobs ~ weibull(alpha, exp(-(mu + Xobs_bg * beta_bg)/alpha));
  increment_log_prob(weibull_ccdf_log(ycen, alpha, exp(-(mu + Xcen_bg * beta_bg)/alpha)));

  beta_bg_raw ~ normal(0.0, 1.0);
  alpha_raw ~ normal(0.0, 1.0);

  mu ~ normal(0.0, tau_mu);
}
