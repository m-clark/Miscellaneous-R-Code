% matlab from Statistical Modeling and Computation (2014 p 314). See the 
% associated twofactorRE.R file for details.

function sfran2_loglike(mu, eta_alpha, eta_gamma, eta, y, Xalpha, Xgamma)
  sigma2_alpha = exp(eta_alpha);
  sigma2_gamma = exp(eta_gamma);
  sigma2 = exp(eta);
  n = length(y);
  Sigma = sigma2*speye(n) + sigma2_alpha*(Xalpha*Xalpha') + sigma2_gamma*(Xgamma*Xgamma');
  l = -n/2*log(2*pi) - sum(log(diag(chol(Sigma)))) - .5*(y-mu)' * (Sigma\(y-mu));
end

y = [1.39 1.29 1.12 1.16 1.52 1.62 1.88 1.87 1.24 1.18 .95 .96 .82 .92 1.18 1.20 1.47 1.41 1.57 1.65];
Xalpha = kron(speye(5), ones(4,1));
Xgamma = kron(speye(10), ones(2,1));

f = @(theta) -sfran_loglike(theta(1), theta(2), theta(3), theta(4), y, Xalpha, Xgamma);
yhat = mean(reshape(y, 4, 5));
theta0 = [mean(y) log(var(yhat)) log(var(y)/3) log(var(y)/3)];
thetahat = fminsearch(f, theta0)
