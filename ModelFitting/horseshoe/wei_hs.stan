/*  Variable naming:
 obs       = observed
 cen       = (right) censored
 N         = number of samples
 M         = number of covariates
 bg        = established risk (or protective) factors
 biom      = candidate biomarkers (candidate risk factors)
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

  vector hs_prior_lp(real r1_global, real r2_global, vector r1_local, vector r2_local) {
    r1_global ~ normal(0.0, 1.0);
    r2_global ~ inv_gamma(0.5, 0.5);

    r1_local ~ normal(0.0, 1.0);
    r2_local ~ inv_gamma(0.5, 0.5);

    return (r1_global * sqrt(r2_global)) * r1_local .* sqrt_vec(r2_local);
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
  int<lower=0> M_biom;
  vector[Nobs] yobs;
  vector[Ncen] ycen;
  matrix[Nobs, M_bg] Xobs_bg;
  matrix[Ncen, M_bg] Xcen_bg;
  matrix[Nobs, M_biom] Xobs_biom;
  matrix[Ncen, M_biom] Xcen_biom;
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

  real<lower=0> tau_s1_biom_raw;
  real<lower=0> tau_s2_biom_raw;
  vector<lower=0>[M_biom] tau1_biom_raw;
  vector<lower=0>[M_biom] tau2_biom_raw;

  real alpha_raw;
  vector[M_bg] beta_bg_raw;
  vector[M_biom] beta_biom_raw;

  real mu;
}

transformed parameters {
  vector[M_biom] beta_biom;
  vector[M_bg] beta_bg;
  real alpha;

  beta_biom <- hs_prior_lp(tau_s1_biom_raw, tau_s2_biom_raw, tau1_biom_raw, tau2_biom_raw) .* beta_biom_raw;
  beta_bg <- bg_prior_lp(tau_s_bg_raw, tau_bg_raw) .* beta_bg_raw;
  alpha <- exp(tau_al * alpha_raw);
}

model {
  yobs ~ weibull(alpha, exp(-(mu + Xobs_bg * beta_bg + Xobs_biom * beta_biom)/alpha));
  increment_log_prob(weibull_ccdf_log(ycen, alpha, exp(-(mu + Xcen_bg * beta_bg + Xcen_biom * beta_biom)/alpha)));

  beta_biom_raw ~ normal(0.0, 1.0);
  beta_bg_raw ~ normal(0.0, 1.0);
  alpha_raw ~ normal(0.0, 1.0);

  mu ~ normal(0.0, tau_mu);
}
